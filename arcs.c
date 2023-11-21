#include "cplex.h"
#include "lapacke.h"
#include "cblas.h"

typedef struct
{
    double start_angle;
    double end_angle;
} Arc;

typedef struct
{
    double x;
    double y;
    double radius;
    Arc *arcs_to_include;
    int included_arc_count;
    bool arc_excluded;
    int max_array_length;
} Circle;

typedef struct
{
    double x_bounds[2];
    double y_bounds[2];
} Box;

enum separation_lp_status { separation_success,
                            separation_failure,
                            separation_nonoptimal,
                            separation_not_strong_enough };

enum exposed_point_status { exposed_points_exist,
                            no_exposed_points_infeasible,
                            no_exposed_points_feasible };

void initialize_circle(double x, double y, double radius, int circle_count, Circle *circle) {
    circle->x = x;
    circle->y = y;
    circle->radius = radius;
    circle->max_array_length = circle_count + 4; // 4 box sides, M other circles
    circle->arcs_to_include = calloc(circle->max_array_length, sizeof(Arc));
    circle->included_arc_count = 1;
    circle->arc_excluded = false;
}

void initialize_box(double x_lb, double x_ub, double y_lb, double y_ub, Box *box) {
    box->x_bounds[0] = x_lb;
    box->x_bounds[1] = x_ub;
    box->y_bounds[0] = y_lb;
    box->y_bounds[1] = y_ub;
}

void free_circle_arcs(Circle *circle) {
    free(circle->arcs_to_include);
}

void _debug_print_arcs(Arc *arcs, int length) {
    if (length < 1) {
        printf("    none\n");
        return;
    }
    for (int i = 0; i < length; i++) {
        printf("    [%.10lf, %.10lf]\n", (arcs + i)->start_angle, (arcs + i)->end_angle);
    }
}

void _debug_print_circle(Circle *circle) {
    printf("Circle with center at (%.10lf, %.10lf), radius of %.10lf\n", circle->x, circle->y, circle->radius);
    printf("    Included arcs:\n");
    if (circle->included_arc_count == 0) {
        printf("        none\n");
        return;
    }
    _debug_print_arcs(circle->arcs_to_include, circle->included_arc_count);
}

void _debug_print_line(int length) {
    for (int i = 0; i < length; i++) {
        printf("=");
    }
    printf("\n");
}

bool _debug_arcs_redundant(Arc *arc1, Arc *arc2) {
    return fabs(arc1->start_angle - arc2->start_angle) < EPSILON && fabs(arc1->end_angle - arc2->end_angle) < EPSILON;
}

bool _debug_arc_already_exists(Arc *arcs, int length, Arc *arc) {
    for (int i = 0; i < length; i++) {
        if (_debug_arcs_redundant(arcs + i, arc)) {
            return true;
        }
    }
    return false;
}

bool _debug_exists_redundant_arc(Arc *arcs, int length) {
    for (int i = 0; i < length; i++) {
        for (int j = i + 1; j < length; j++) {
            if (_debug_arcs_redundant(arcs + i, arcs + j)) {
                return true;
            }
        }
    }
    return false;
}

bool determine_circle_exclusion(Box *box, double x, double y, double radius) {
    /*
    Returns true if we determine that the circle centered at (x,y) is entirely outside the box
    */
    return (x + radius < box->x_bounds[0]) || (x - radius > box->x_bounds[1]) || (y + radius < box->y_bounds[0]) || (y - radius > box->y_bounds[1]);
}

void find_k_smallest(double *arr, int length, int k, double *k_smallest, int *indices_smallest) {
    /*
    Writes k smallest elements in increasing order
    */
    for (int i = 0; i < k; i++) {
        k_smallest[i] = INFINITY;
        indices_smallest[i] = -1;
    }
    for (int i = 0; i < length; i++) {
        int insertion_index = 0;
        while (arr[i] > k_smallest[insertion_index]) {
            insertion_index++;
            if (insertion_index == k) {
                break;
            }
        }
        if (insertion_index == k) {
            continue;
        }
        if (insertion_index == k - 1) {
            k_smallest[insertion_index] = arr[i];
            indices_smallest[insertion_index] = i;
            continue;
        }
        double temp = k_smallest[insertion_index];
        int temp_i = indices_smallest[insertion_index];
        k_smallest[insertion_index] = arr[i];
        indices_smallest[insertion_index] = i;
        while (insertion_index < k - 1) {
            double temp2 = k_smallest[insertion_index + 1];
            int temp_i2 = indices_smallest[insertion_index + 1];
            k_smallest[insertion_index + 1] = temp;
            indices_smallest[insertion_index + 1] = temp_i;
            temp = temp2;
            temp_i = temp_i2;
            insertion_index++;
        }
    }
}

void copy_arc(Arc *source, Arc *destination) {
    destination->start_angle = source->start_angle;
    destination->end_angle = source->end_angle;
}

double calculate_angle(double x, double y) {
    /*
    Angle of (x, y) with respect to the origin in radians
    */
    double arctan = atan2(y, x);
    if (y < 0) {
        arctan += TWO_PI;
    }
    return arctan;
}

double solve_quadratic(double a, double b, double c, double x_lb, double x_ub) {
    /*
    Assumes nonnegative discriminant. Returns solution in the specified interval
    */
    double discriminant_scaled = sqrt(pow(b, 2) - 4 * a * c) / (2 * a);
    double answer = (-b) / (2 * a) + discriminant_scaled;
    if (answer >= x_lb && answer <= x_ub) {
        return answer;
    }
    return answer - 2 * discriminant_scaled;
}

void two_root_quadratic(double a, double b, double c, double *root1, double *root2) {
    /*
    Assumes nonnegative discriminant. Calculates both solutions
    */
    double second_term_scaled = sqrt(pow(b, 2) - 4 * a * c) / (2 * a);
    *(root1) = (-b) / (2 * a) + second_term_scaled;
    *(root2) = *(root1)-2 * second_term_scaled;
}

void points_to_arc(double x1, double y1, double x2, double y2, double reference_x, double reference_y, Arc *arc) {
    /*
    Given two points (x1, y1) and (x2, y2) on a shape's boundary,
    calculates the starting and ending angles of the boundary segment defined by starting at (x1, y1)
    and going counterclockwise to (x2, y2).
    Must be provided with a reference point to base the angle calculations on.

    Angles are always less than 2*pi, so if, for example,
    the starting angle is 5*pi/6 and the ending angle is pi/6, we change the starting angle
    to be -pi/6

    Also stores the index of the other shape in the intersection
    */
    double theta1 = calculate_angle(x1 - reference_x, y1 - reference_y);
    double theta2 = calculate_angle(x2 - reference_x, y2 - reference_y);
    if (theta1 > theta2) {
        theta1 -= TWO_PI;
    }
    arc->start_angle = theta1;
    arc->end_angle = theta2;
}

bool arc_within_arc(Arc *a, Arc *b) {
    /*
    Checks whether a \subseteq b assuming they belong to the same circle
    */
    double a_start_minus2pi = a->start_angle - TWO_PI;
    double a_end_minus2pi = a->end_angle - TWO_PI;
    return (
        ((a->start_angle >= b->start_angle) && (a->end_angle <= b->end_angle)) ||
        ((a_start_minus2pi >= b->start_angle) && (a_end_minus2pi <= b->end_angle)));
}

bool angle_within_arc(double angle, Arc *arc) {
    double angle_minus2pi = angle - TWO_PI;
    double angle_plus2pi = angle + TWO_PI;
    return (
        (angle >= arc->start_angle && angle <= arc->end_angle) ||
        (angle_minus2pi >= arc->start_angle && angle_minus2pi <= arc->end_angle) ||
        (angle_plus2pi >= arc->start_angle && angle_plus2pi <= arc->end_angle));
}

void standardize_arc_angles(Arc *arc) {
    /*
    Assumes the points at the start and end angles are correct, but a factor of 2*pi may be incorrect.
    In particular, we can assume the given angles are all in [-2*pi, 2*pi]
    Changes the angles so that the arc length is <= 2*pi, the start angle is in [-2*pi, 2*pi], the
    end angle is in [0, 2*pi], and the end angle is greater than the start angle
    */
    if (arc->end_angle < 0) {
        arc->end_angle += TWO_PI;
    }
    if (arc->start_angle > arc->end_angle) {
        arc->start_angle -= TWO_PI;
    }
    if (arc->end_angle - arc->start_angle > TWO_PI) {
        arc->start_angle += TWO_PI;
    }
}

void update_included_arcs(Arc *arcs_to_include, int *included_arc_count, bool *arc_excluded, Arc *arc_to_exclude) {
    /*
    Calculates effect of excluding the given arc on the included arcs.
    */
    // if no arcs have been excluded yet, the single included arc is the complement of the given arc
    if (!(*arc_excluded)) {
        /*
        arcs to exclude are never [0, 2*pi], so we don't need to check that case here
        we are also guaranteed that all angles are within [-2*pi, 2*pi], start angles are only negative if
        0 is in the arc, and end angles are always nonnegative
        */
        if (arc_to_exclude->start_angle < 0) {
            arcs_to_include[0].start_angle = arc_to_exclude->end_angle;
            arcs_to_include[0].end_angle = arc_to_exclude->start_angle + TWO_PI;
        } else {
            arcs_to_include[0].start_angle = arc_to_exclude->end_angle - TWO_PI;
            arcs_to_include[0].end_angle = arc_to_exclude->start_angle;
        }
        *arc_excluded = true;
        return;
    }
    /*
    The arc to exclude might contain several included arcs. We mark their indices and then copy over the
    non-deleted (but possibly modified) arcs at the end. Note that this is mutually exclusive with the arc splitting.
    */
    bool deleted[*included_arc_count];
    for (int i = 0; i < *included_arc_count; i++) {
        // case 1: excluded arc entirely within an included arc => split that arc and we're done
        if (arc_within_arc(arc_to_exclude, arcs_to_include + i)) {
            // split and append the new arc to the end
            double old_end_angle = arcs_to_include[i].end_angle;
            if (arc_to_exclude->start_angle < 0) {
                arcs_to_include[i].start_angle += TWO_PI;
                arcs_to_include[i].end_angle = arc_to_exclude->start_angle + TWO_PI;
            } else {
                arcs_to_include[i].end_angle = arc_to_exclude->start_angle;
            }
            arcs_to_include[*included_arc_count].start_angle = arc_to_exclude->end_angle;
            arcs_to_include[*included_arc_count].end_angle = old_end_angle;
            *included_arc_count += 1;
            return;
        }
        // case 2: excluded arc covers this arc => delete this arc. This might be true for multiple arcs
        if (arc_within_arc(arcs_to_include + i, arc_to_exclude)) {
            deleted[i] = true;
            continue;
        }
        deleted[i] = false;
        bool angles_may_need_standardizing = false;
        // case 3: an endpoint of this arc is within the excluded arc
        if (angle_within_arc(arcs_to_include[i].start_angle, arc_to_exclude)) {
            arcs_to_include[i].start_angle = arc_to_exclude->end_angle;
            angles_may_need_standardizing = true;
        }
        if (angle_within_arc(arcs_to_include[i].end_angle, arc_to_exclude)) {
            arcs_to_include[i].end_angle = arc_to_exclude->start_angle;
            angles_may_need_standardizing = true;
        }
        if (angles_may_need_standardizing) {
            standardize_arc_angles(arcs_to_include + i);
        }
    }
    int arcs_copied = 0;
    for (int i = 0; i < *included_arc_count; i++) {
        if (deleted[i]) {
            continue;
        }
        copy_arc(arcs_to_include + i, arcs_to_include + arcs_copied);
        arcs_copied++;
    }
    *included_arc_count = arcs_copied;
}

void calculate_box_intersections(Circle *circles, int circle_index, int circle_count, Box *box) {
    /*
    Updates included arcs based on regions of each circle boundary that are outside the box
    */
    Circle *circle = circles + circle_index;
    if (circle->included_arc_count == 0) {
        return;
    }
    double x_lb = box->x_bounds[0];
    double x_ub = box->x_bounds[1];
    double y_lb = box->y_bounds[0];
    double y_ub = box->y_bounds[1];
    double x1, y1, x2, y2, linear_coefficient, constant_term;
    Arc *arc = calloc(1, sizeof(Arc));
    if (circle->x + circle->radius > x_ub) {
        // point 1 in [y_lb, circle y], point 2 with y in [circle y, y_ub]
        linear_coefficient = -2 * circle->y;
        constant_term = pow(circle->y, 2) + pow(x_ub - circle->x, 2) - pow(circle->radius, 2);
        x1 = x_ub;
        y1 = solve_quadratic(1, linear_coefficient, constant_term, y_lb, circle->y);
        x2 = x_ub;
        y2 = solve_quadratic(1, linear_coefficient, constant_term, circle->y, y_ub);
        points_to_arc(x1, y1, x2, y2, circle->x, circle->y, arc);
        update_included_arcs(circle->arcs_to_include, &(circle->included_arc_count), &(circle->arc_excluded), arc);
    }
    if (circle->x - circle->radius < x_lb) {
        // point 1 with y in [circle y, y_ub], point 2 in [y_lb, circle y]
        linear_coefficient = -2 * circle->y;
        constant_term = pow(circle->y, 2) + pow(x_lb - circle->x, 2) - pow(circle->radius, 2);
        x1 = x_lb;
        y1 = solve_quadratic(1, linear_coefficient, constant_term, circle->y, y_ub);
        x2 = x_lb;
        y2 = solve_quadratic(1, linear_coefficient, constant_term, y_lb, circle->y);
        points_to_arc(x1, y1, x2, y2, circle->x, circle->y, arc);
        update_included_arcs(circle->arcs_to_include, &(circle->included_arc_count), &(circle->arc_excluded), arc);
    }
    if (circle->y + circle->radius > y_ub) {
        // point 1 with x in [circle x, x_ub], point 2 in [x_lb, circle x]
        linear_coefficient = -2 * circle->x;
        constant_term = pow(circle->x, 2) + pow(y_ub - circle->y, 2) - pow(circle->radius, 2);
        x1 = solve_quadratic(1, linear_coefficient, constant_term, circle->x, x_ub);
        y1 = y_ub;
        x2 = solve_quadratic(1, linear_coefficient, constant_term, x_lb, circle->x);
        y2 = y_ub;
        points_to_arc(x1, y1, x2, y2, circle->x, circle->y, arc);
        update_included_arcs(circle->arcs_to_include, &(circle->included_arc_count), &(circle->arc_excluded), arc);
    }
    if (circle->y - circle->radius < y_lb) {
        // point 1 in [x_lb, circle x], point 2 with x in [circle x, x_ub]
        linear_coefficient = -2 * circle->x;
        constant_term = pow(circle->x, 2) + pow(y_lb - circle->y, 2) - pow(circle->radius, 2);
        x1 = solve_quadratic(1, linear_coefficient, constant_term, x_lb, circle->x);
        y1 = y_lb;
        x2 = solve_quadratic(1, linear_coefficient, constant_term, circle->x, x_ub);
        y2 = y_lb;
        points_to_arc(x1, y1, x2, y2, circle->x, circle->y, arc);
        update_included_arcs(circle->arcs_to_include, &(circle->included_arc_count), &(circle->arc_excluded), arc);
    }
    free(arc);
}

void calculate_circle_intersection(Circle *circles, int index1, int index2) {
    /*
    If two circle boundaries intersect at two distinct points, we add the corresponding arcs to the
    list of arcs to exclude for each circle
    */
    Circle *circle1 = circles + index1;
    Circle *circle2 = circles + index2;
    if (pow(circle1->x - circle2->x, 2) + pow(circle1->y - circle2->y, 2) >= pow(circle1->radius + circle2->radius, 2)) {
        // intersection of circle boundaries is empty or is a singleton
        return;
    }
    if (circle1->included_arc_count == 0 || circle2->included_arc_count == 0) {
        // every point on one circle's boundary is covered by another cirle, making this calculation redundant
        return;
    }
    double linear_coefficient, constant_term;
    Arc arc[1];
    // if the line through the circle centers is vertical, we avoid calculating its slope
    if (fabs(circle1->x - circle2->x) < EPSILON) {
        double y = (pow(circle1->radius, 2) - pow(circle2->radius, 2) +
                    pow(circle2->y, 2) - pow(circle1->y, 2)) /
                   (-2 * (circle1->y - circle2->y));
        double x1, x2;
        linear_coefficient = -2 * circle1->x;
        constant_term = pow(circle1->x, 2) + pow(y - circle1->y, 2) - pow(circle1->radius, 2);
        two_root_quadratic(1, linear_coefficient, constant_term, &x1, &x2);
        // assign arcs and add to circles
        if (circle1->y < circle2->y) {
            if (x1 < x2) {
                points_to_arc(x2, y, x1, y, circle1->x, circle1->y, arc);
                update_included_arcs(circle1->arcs_to_include, &(circle1->included_arc_count), &(circle1->arc_excluded), arc);
                points_to_arc(x1, y, x2, y, circle2->x, circle2->y, arc);
                update_included_arcs(circle2->arcs_to_include, &(circle2->included_arc_count), &(circle2->arc_excluded), arc);
            } else {
                points_to_arc(x1, y, x2, y, circle1->x, circle1->y, arc);
                update_included_arcs(circle1->arcs_to_include, &(circle1->included_arc_count), &(circle1->arc_excluded), arc);
                points_to_arc(x2, y, x1, y, circle2->x, circle2->y, arc);
                update_included_arcs(circle2->arcs_to_include, &(circle2->included_arc_count), &(circle2->arc_excluded), arc);
            }
        } else {
            if (x1 < x2) {
                points_to_arc(x2, y, x1, y, circle2->x, circle2->y, arc);
                update_included_arcs(circle2->arcs_to_include, &(circle2->included_arc_count), &(circle2->arc_excluded), arc);
                points_to_arc(x1, y, x2, y, circle1->x, circle1->y, arc);
                update_included_arcs(circle1->arcs_to_include, &(circle1->included_arc_count), &(circle1->arc_excluded), arc);
            } else {
                points_to_arc(x1, y, x2, y, circle2->x, circle2->y, arc);
                update_included_arcs(circle2->arcs_to_include, &(circle2->included_arc_count), &(circle2->arc_excluded), arc);
                points_to_arc(x2, y, x1, y, circle1->x, circle1->y, arc);
                update_included_arcs(circle1->arcs_to_include, &(circle1->included_arc_count), &(circle1->arc_excluded), arc);
            }
        }
        return;
    }
    /*
    For circles centered at (a1,b1) and (a2,b2), we solve the system of equations
    (x - a1)^2 + (y - b1)^2 = r1^2
    (x - a2)^2 + (y - b2)^2 = r2^2
    Subtracting the second from the first, the quadratic terms disappear. After rearranging, we get
    y = alpha*x - beta, where
    alpha = (a2 - a1) / (b1 - b2), and
    beta = -(r1^2 - r2^2 + a2^2 - a1^2 + b2^2 - b1^2) / (2 * (b1 - b2))
    We can substitute this back into the first quadratic equation to get
    p0 x^2 + p1 x + p2 = 0, where
    p0 = 1 + alpha^2
    p1 = -2 * (a1 + alpha*beta - b1*alpha)
    p2 = a1^2 + b1^2 + beta^2 + 2*b1*beta - r1^2
    after which we apply the quadratic formula.
    */
    double alpha = (circle2->x - circle1->x) / (circle1->y - circle2->y);
    double beta = (pow(circle1->radius, 2) - pow(circle2->radius, 2) +
                   pow(circle2->x, 2) - pow(circle1->x, 2) +
                   pow(circle2->y, 2) - pow(circle1->y, 2)) /
                  (2 * (circle1->y - circle2->y));
    double p0 = 1 + pow(alpha, 2);
    double p1 = -2 * (circle1->x + alpha * beta + (circle1->y) * alpha);
    double p2 = pow(circle1->x, 2) + pow(circle1->y, 2) + pow(beta, 2) + 2 * (circle1->y) * beta - pow(circle1->radius, 2);
    double discriminant = pow(p1, 2) - 4 * p0 * p2;
    if (discriminant < 0) {
        return;
    }
    double second_term_scaled = sqrt(pow(p1, 2) - 4 * p0 * p2) / (2 * p0);
    double first_term = -p1 / (2 * p0);
    double x1 = first_term + second_term_scaled;
    double y1 = alpha * x1 - beta;
    double x2 = first_term - second_term_scaled;
    double y2 = alpha * x2 - beta;

    points_to_arc(x1, y1, x2, y2, circle1->x, circle1->y, arc);
    double theta_21 = calculate_angle(circle2->x - circle1->x, circle2->y - circle1->y);
    if (angle_within_arc(theta_21, arc)) {
        update_included_arcs(circle1->arcs_to_include, &(circle1->included_arc_count), &(circle1->arc_excluded), arc);
        points_to_arc(x2, y2, x1, y1, circle2->x, circle2->y, arc);
        update_included_arcs(circle2->arcs_to_include, &(circle2->included_arc_count), &(circle2->arc_excluded), arc);
        return;
    }
    // reverse start and end angles
    double new_start_angle;
    double new_end_angle;
    if (arc->start_angle >= 0) {
        new_start_angle = arc->end_angle - TWO_PI;
        new_end_angle = arc->start_angle;
    } else {
        new_start_angle = arc->end_angle;
        new_end_angle = arc->start_angle + TWO_PI;
    }
    arc->start_angle = new_start_angle;
    arc->end_angle = new_end_angle;
    update_included_arcs(circle1->arcs_to_include, &(circle1->included_arc_count), &(circle1->arc_excluded), arc);
    points_to_arc(x1, y1, x2, y2, circle2->x, circle2->y, arc);
    update_included_arcs(circle2->arcs_to_include, &(circle2->included_arc_count), &(circle2->arc_excluded), arc);
}

void angle_to_point(Circle *circle, double theta, double *x, double *y) {
    double x_offset = circle->radius * cos(theta);
    *(x) = circle->x + x_offset;
    double y_offset = circle->radius * sin(theta);
    *(y) = circle->y + y_offset;
}

bool point_in_array(double *x_arr, double *y_arr, int length, double x, double y) {
    for (int i = 0; i < length; i++) {
        if (fabs(x - *(x_arr + i)) + fabs(y - *(y_arr + i)) < EPSILON) {
            return true;
        }
    }
    return false;
}

bool point_is_valid(double *x_exposed, double *y_exposed, int exposed_point_count,
                    double x_lb, double x_ub, double y_lb, double y_ub,
                    double x, double y, Circle *circles, int circle_count) {
    /*
    Checks whether point is already in the collection of exposed points, outside the box,
    or within another circle
    */
    if (x < x_lb - EPSILON || x > x_ub + EPSILON || y < y_lb - EPSILON || y > y_ub + EPSILON) {
        return false;
    }
    if (point_in_array(x_exposed, y_exposed, exposed_point_count, x, y)) {
        return false;
    }
    for (int i = 0; i < circle_count; i++) {
        if (pow(circles[i].x - x, 2) + pow(circles[i].y - y, 2) < pow(circles[i].radius, 2) - EPSILON) {
            return false;
        }
    }
    return true;
}

void write_arc_endpoints(Circle *circles, int circle_count, int circle_index, double x_lb, double x_ub, double y_lb, double y_ub,
                         Arc *arcs_to_include, int included_arc_count,
                         double *x_exposed, double *y_exposed, int *exposed_point_count) {
    /*
    Translates arc endpoints to Cartesian coordinates and checks for uniqueness and correctness
    */
    int added_point_count = 0;
    Circle *circle = circles + circle_index;
    // if no arcs were excluded, exit here rather than using the included arc [0, 2*pi]
    if (!circle->arc_excluded) {
        return;
    }
    for (int i = 0; i < included_arc_count; i++) {
        double x;
        double y;
        angle_to_point(circle, (arcs_to_include + i)->start_angle, &x, &y);
        if (point_is_valid(x_exposed, y_exposed, *(exposed_point_count), x_lb, x_ub, y_lb, y_ub, x, y, circles, circle_count)) {
            *(x_exposed + *(exposed_point_count)) = x;
            *(y_exposed + *(exposed_point_count)) = y;
            *(exposed_point_count) += 1;
            added_point_count++;
        }
        angle_to_point(circle, (arcs_to_include + i)->end_angle, &x, &y);
        if (point_is_valid(x_exposed, y_exposed, *(exposed_point_count), x_lb, x_ub, y_lb, y_ub, x, y, circles, circle_count)) {
            *(x_exposed + *(exposed_point_count)) = x;
            *(y_exposed + *(exposed_point_count)) = y;
            *(exposed_point_count) += 1;
            added_point_count++;
        }
    }
}

int collect_exposed_points(Circle circles[], int circle_count, double x_lb, double x_ub,
                           double y_lb, double y_ub, double *x_exposed, double *y_exposed, int *max_exposed_points) {
    /*
    Calculates unique exposed points of the region defined by the box [x_lb, x_ub] x [y_lb, y_ub] without the circles,
    given the arcs to exclude for each circle
    Changes the number of max exposed points to be the actual number of exposed points.
    */
    int exposed_point_count = 0;
    int total_included_arcs = 0;
    for (int i = 0; i < circle_count; i++) {
        if (circles[i].included_arc_count == 0) {
            continue;
        }
        total_included_arcs += circles[i].included_arc_count;
        write_arc_endpoints(circles, circle_count, i, x_lb, x_ub, y_lb, y_ub, circles[i].arcs_to_include, circles[i].included_arc_count, x_exposed, y_exposed, &exposed_point_count);
    }
    if (total_included_arcs == 0) {
        return no_exposed_points_infeasible;
    }
    // check box corners
    double x_bounds[2] = {x_lb, x_ub};
    double y_bounds[2] = {y_lb, y_ub};
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            double x_bound = x_bounds[i];
            double y_bound = y_bounds[j];
            bool corner_in_circle = false;
            for (int k = 0; k < circle_count; k++) {
                if (pow(x_bound - circles[k].x, 2) + pow(y_bound - circles[k].y, 2) <= pow(circles[k].radius, 2)) {
                    corner_in_circle = true;
                    break;
                }
            }
            if (!corner_in_circle) {
                *(x_exposed + exposed_point_count) = x_bound;
                *(y_exposed + exposed_point_count) = y_bound;
                exposed_point_count++;
            }
        }
    }
    *(max_exposed_points) = exposed_point_count;
    if (exposed_point_count == 0) {
        return no_exposed_points_feasible;
    }
    return exposed_points_exist;
}

void print_cplex_error(CPXENVptr env, int cplex_status) {
    char errmsg[CPXMESSAGEBUFSIZE];
    CPXgeterrorstring(env, cplex_status, errmsg);
    fprintf(stderr, "%s", errmsg);
}

int initialize_lp(CPXENVptr *env, CPXLPptr *lp) {
    int cplex_status = 0;
    *(env) = CPXopenCPLEX(&cplex_status);
    if (*(env) == NULL) {
        fprintf(stderr, "Could not open CPLEX environment.\n");
        print_cplex_error(*(env), cplex_status);
        return 1;
    }
    cplex_status = CPXsetintparam(*(env), CPXPARAM_Read_DataCheck, CPX_DATACHECK_ASSIST);
    if (cplex_status) {
        fprintf(stderr, "Failure to turn on data checking, error %d.\n", cplex_status);
        return 1;
    }
    *(lp) = CPXcreateprob(*(env), &cplex_status, "separation");
    if (*(lp) == NULL) {
        fprintf(stderr, "Failed to create LP.\n");
        return 1;
    }
    cplex_status = CPXchgobjsen(*(env), *(lp), CPX_MAX);
    if (cplex_status) {
        fprintf(stderr, "Failed to set direction to maximization.\n");
        return 1;
    }
    return 0;
}

void clean_lp(CPXENVptr *env, CPXLPptr *lp) {
    /*
    If we are unable to free the CPLEX environment, we should halt execution altogether to avoid
    memory leaks and crashes.
    */
    int cplex_status = 0;
    if (*(lp) == NULL) {
        fprintf(stderr, "Fatal error: attempted to free null CPLEX LP.\n");
        exit(1);
    }
    cplex_status = CPXfreeprob(*(env), lp);
    if (cplex_status) {
        fprintf(stderr, "Fatal error: CPXfreeprob failed, error code %d.\n", cplex_status);
        exit(1);
    }
    if (*(env) == NULL) {
        fprintf(stderr, "Fatal error: attempted to free null CPLEX environment.\n");
        exit(1);
    }
    cplex_status = CPXcloseCPLEX(env);
    if (cplex_status) {
        fprintf(stderr, "Fatal error: Could not close CPLEX environment.\n");
        print_cplex_error(*(env), cplex_status);
        exit(1);
    }
}

int solve_separation_lp(double *x_exposed, double *y_exposed, int exposed_point_count, double x, double y, double infinity,
                        double *x_coefficient, double *y_coefficient, double *constant_term) {
    /*
    Separate (x,y) from the exposed points, if possible
    */
    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    int initialization_status = initialize_lp(&env, &lp);
    if (initialization_status != 0) {
        return separation_failure;
    }
    int cplex_status = 0;
    double objective_coefficients[3] = {x, y, 1};
    double lower_bounds[3] = {-1.0, -1.0, -1.0};
    double upper_bounds[3] = {1.0, 1.0, 1.0};
    char *column_names[3] = {"x_coefficient", "y_coefficient", "constant_term"};
    cplex_status = CPXnewcols(env, lp, 3, objective_coefficients, lower_bounds, upper_bounds, NULL, column_names);
    if (cplex_status) {
        fprintf(stderr, "Could not add variables to CPLEX LP.\n");
        print_cplex_error(env, cplex_status);
        clean_lp(&env, &lp);
        return separation_failure;
    }
    int nonzero_count = 3 * exposed_point_count;
    int row_starts[exposed_point_count];
    int columns[nonzero_count];
    double coefficients[nonzero_count];
    char constraint_sense[exposed_point_count];
    double right_hand_sides[exposed_point_count];
    for (int i = 0; i < exposed_point_count; i++) {
        row_starts[i] = 3 * i;
        coefficients[3 * i] = *(x_exposed + i);
        columns[3 * i] = 0;
        coefficients[3 * i + 1] = *(y_exposed + i);
        columns[3 * i + 1] = 1;
        coefficients[3 * i + 2] = 1.0;
        columns[3 * i + 2] = 2;
        constraint_sense[i] = 'L';
        right_hand_sides[i] = 0.0;
    }
    cplex_status = CPXaddrows(env, lp, 0, exposed_point_count, nonzero_count, right_hand_sides, constraint_sense,
                              row_starts, columns, coefficients, NULL, NULL);
    if (cplex_status) {
        fprintf(stderr, "Could not add constraints to CPLEX LP.\n");
        print_cplex_error(env, cplex_status);
        clean_lp(&env, &lp);
        return separation_failure;
    }
    cplex_status = CPXlpopt(env, lp);
    if (cplex_status) {
        fprintf(stderr, "Failed to solve LP.\n");
        print_cplex_error(env, cplex_status);
        clean_lp(&env, &lp);
        return separation_failure;
    }
    int solution_status;
    double solution[3] = {0, 0, 0};
    double objective_value;
    cplex_status = CPXsolution(env, lp, &solution_status, &objective_value, solution, NULL, NULL, NULL);
    if (cplex_status) {
        clean_lp(&env, &lp);
        return separation_failure;
    }
    if (solution_status != CPX_STAT_OPTIMAL) {
        fprintf(stderr, "LP not solved to optimality; code %d\n", solution_status);
        clean_lp(&env, &lp);
        return separation_nonoptimal;
    }
    double cut_violation_tolerance = 1e-4;
    if (objective_value < cut_violation_tolerance) {
        clean_lp(&env, &lp);
        return separation_not_strong_enough;
    }

    *(x_coefficient) = solution[0];
    *(y_coefficient) = solution[1];
    *(constant_term) = -solution[2];
    clean_lp(&env, &lp);
    return separation_success;
}

int solve_concave_envelope_lp(double *x_exposed, double *y_exposed, double *sq_distances, int exposed_point_count,
                              double x, double y, double d, double infinity, double *optimal_coefficients) {
    /*
    Generate a facet of the concave envelope of the binding distance, if possible
    */
    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    int initialization_status = initialize_lp(&env, &lp);
    if (initialization_status != 0) {
        return separation_failure;
    }
    int cplex_status = 0;
    double objective_coefficients[exposed_point_count];
    double lower_bounds[exposed_point_count];
    double upper_bounds[exposed_point_count];
    for (int i = 0; i < exposed_point_count; i++) {
        objective_coefficients[i] = sq_distances[i];
        lower_bounds[i] = 0;
        upper_bounds[i] = 1;
    }
    cplex_status = CPXnewcols(env, lp, exposed_point_count, objective_coefficients, lower_bounds, upper_bounds, NULL, NULL);
    if (cplex_status) {
        fprintf(stderr, "Could not add variables to CPLEX LP.\n");
        print_cplex_error(env, cplex_status);
        clean_lp(&env, &lp);
        return separation_failure;
    }
    int nonzero_count = 3 * exposed_point_count;
    int row_starts[3] = {0, exposed_point_count, 2 * exposed_point_count};
    int columns[nonzero_count];
    double coefficients[nonzero_count];
    char constraint_senses[3] = {'E', 'E', 'E'};
    double right_hand_sides[3] = {x, y, 1.0};
    for (int i = 0; i < exposed_point_count; i++) {
        coefficients[i] = *(x_exposed + i);
        columns[i] = i;
        coefficients[exposed_point_count + i] = *(y_exposed + i);
        columns[exposed_point_count + i] = i;
        coefficients[2 * exposed_point_count + i] = 1.0;
        columns[2 * exposed_point_count + i] = i;
    }
    cplex_status = CPXaddrows(env, lp, 0, 3, nonzero_count, right_hand_sides, constraint_senses,
                              row_starts, columns, coefficients, NULL, NULL);
    if (cplex_status) {
        fprintf(stderr, "Could not add constraints to CPLEX LP.\n");
        print_cplex_error(env, cplex_status);
        clean_lp(&env, &lp);
        return separation_failure;
    }
    cplex_status = CPXlpopt(env, lp);
    if (cplex_status) {
        fprintf(stderr, "Failed to solve LP.\n");
        print_cplex_error(env, cplex_status);
        clean_lp(&env, &lp);
        return separation_failure;
    }
    double objective_value;
    int solution_status;
    cplex_status = CPXsolution(env, lp, &solution_status, &objective_value, optimal_coefficients, NULL, NULL, NULL);
    if (cplex_status) {
        clean_lp(&env, &lp);
        return separation_failure;
    }
    if (solution_status != CPX_STAT_OPTIMAL && solution_status != CPX_STAT_INFEASIBLE) {
        fprintf(stderr, "LP not solved to optimality; code %d\n", solution_status);
        clean_lp(&env, &lp);
        return separation_nonoptimal;
    }
    double cut_violation_tolerance = 1e-4;
    if (d < objective_value + cut_violation_tolerance) {
        clean_lp(&env, &lp);
        return separation_not_strong_enough;
    }
    clean_lp(&env, &lp);
    return separation_success;
}

int calculate_plane_coefficients(double *x, double *y, double *z, double *coefficients, double *constant_term) {
    /*
    Given x, y, and z as 3-element arrays, calculates coefficients (a, b, c) of the plane passing through 
    the three points (x_i, y_i, z_i) with equation a*x_i + b*y_i + c*z_i = 1.
    */
    lapack_int n = 3;
    double A[9];
    double b[3] = {1.0, 1.0, 1.0};
    int row_permutation[3];
    for (int i = 0; i < 3; i++) {
        A[3 * i] = x[i];
        A[3 * i + 1] = y[i];
        A[3 * i + 2] = z[i];
    }
    lapack_int status = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, A, n, row_permutation, b, 1);
    if (status != 0) {
        return status;
    }
    *constant_term = 1.0;
    for (int i = 0; i < 3; i++) {
        coefficients[i] = b[i];
    }
    return 0;
}

int calculate_coefficients_from_multipliers(double *x_exposed, double *y_exposed, double *sq_distances, double *optimal_multipliers,
    int point_count, double *cut_coefficients, double *cut_right_hand_side) {
    /*
    Using solution to the concave envelope representation problem, calculates supporting hyperplane to the concave envelope
    */
    double active_x[3];
    double active_y[3];
    double active_d[3];
    int points_collected = 0;
    int i = 0;
    while (points_collected < 3) {
        if (i >= point_count) {
            return -1;
        }
        if (optimal_multipliers[i] > EPSILON) {
            active_x[points_collected] = x_exposed[i];
            active_y[points_collected] = y_exposed[i];
            active_d[points_collected] = sq_distances[i];
            points_collected++;
        }
        i++;
    }
    int status = calculate_plane_coefficients(active_x, active_y, active_d, cut_coefficients, cut_right_hand_side);
    if (status != 0) {
        return status;
    }
    // cuts will always have a positive coefficient for d if written in LHS <= RHS form
    if (cut_coefficients[2] > 0) {
        return 0;
    }
    for (int j = 0; j < 3; j++) {
        cut_coefficients[j] = -cut_coefficients[j];
    }
    *(cut_right_hand_side) = - *(cut_right_hand_side);
    return 0;
}

void free_exposed_points_and_circles(double *x_exposed, double *y_exposed, Circle *circles, int circle_count) {
    free(x_exposed);
    free(y_exposed);
    for (int i = 0; i < circle_count; i++) {
        free_circle_arcs(circles + i);
    }
}