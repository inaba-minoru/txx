#ifndef BEZIER_H
#define BEZIER_H

#include <vecmath.h>

#include <material.hpp>
#include <object3d.hpp>
#include <vector>

#define eps 1e-6
#define INF (1 << 20)
#define sqr(x) ((x) * (x))

class BezierCurve {
   public:
    double *dx, *dy, max, height, max2, r, num;

    int n;

    struct D {
        double t0, t1, width, y0, y1, width2;
    } data[20];

    BezierCurve(std::vector<Vector3f>& p) {
        n = p.size();
        num = n;
        r = 0.365;

        dx = new double[n];
        dy = new double[n];
        assert(std::abs(p[0].y()) <= 1e-6);
        --n;
        // preproces
        for (int i = 0; i <= n; ++i) {
            dx[i] = p[0].x();
            dy[i] = p[0].y();
            for (int j = 0; j <= n - i; ++j) {
                p[j].x() = p[j + 1].x() - p[j].x();
                p[j].y() = p[j + 1].y() - p[j].y();
            }
        }
        double n_down = 1, fac = 1, nxt = n;
        for (int i = 0; i <= n; ++i, --nxt) {
            fac = fac * (i == 0 ? 1 : i);
            dx[i] = dx[i] * n_down / fac;
            dy[i] = dy[i] * n_down / fac;
            n_down *= nxt;
        }
        max = 0;
        double interval = 1. / (num - 1), c = 0;
        for (int cnt = 0; cnt <= num; c += interval, ++cnt) {
            data[cnt].width = 0;
            data[cnt].t0 = std::max(0., c - r);
            data[cnt].t1 = std::min(1., c + r);
            data[cnt].y0 = getpos(data[cnt].t0).y();
            data[cnt].y1 = getpos(data[cnt].t1).y();
            for (double t = data[cnt].t0; t <= data[cnt].t1; t += 0.00001) {
                Vector3f pos = getpos(t);
                if (data[cnt].width < pos.x()) data[cnt].width = pos.x();
            }
            if (max < data[cnt].width) max = data[cnt].width;
            data[cnt].width += eps;
            data[cnt].width2 = sqr(data[cnt].width);
        }
        max += eps;
        max2 = max * max;
        height = getpos(1).y();
    }
    Vector3f getpos(double t) {
        double ans_x = 0, ans_y = 0, t_pow = 1;
        for (int i = 0; i <= n; ++i) {
            ans_x += dx[i] * t_pow;
            ans_y += dy[i] * t_pow;
            t_pow *= t;
        }
        return Vector3f(ans_x, ans_y, 0);
    }
    Vector3f getdir(double t) {
        double ans_x = 0, ans_y = 0, t_pow = 1;
        for (int i = 1; i <= n; ++i) {
            ans_x += dx[i] * i * t_pow;
            ans_y += dy[i] * i * t_pow;
            t_pow *= t;
        }
        return Vector3f(ans_x, ans_y, 0);
    }
    Vector3f getdir2(double t) {
        double ans_x = 0, ans_y = 0, t_pow = 1;
        for (int i = 2; i <= n; ++i) {
            ans_x += dx[i] * i * (i - 1) * t_pow;
            ans_y += dy[i] * i * (i - 1) * t_pow;
            t_pow *= t;
        }
        return Vector3f(ans_x, ans_y, 0);
    }
};

class BezierRotator : public Object3D {
   public:
    BezierCurve curve;
    Vector3f pos;  // the buttom center point
    BezierRotator(Vector3f pos_, BezierCurve c_, Material* m)
        : pos(pos_), curve(c_), Object3D(m) {}
    double solve_t(double yc) {  // solve y(t)=yc
        // assert(0 <= yc && yc <= curve.height);
        double t = .5, ft, dft;
        for (int i = 10; i--;) {
            if (t < 0)
                t = 0;
            else if (t > 1)
                t = 1;
            ft = curve.getpos(t).y() - yc, dft = curve.getdir(t).y();
            if (std::abs(ft) < eps) return t;
            t -= ft / dft;
        }
        return -1;
    }
    virtual Vector3f change_for_bezier(Vector3f inter_p) {
        double t = solve_t(inter_p.y() - pos.y());
        double u = atan2(inter_p.z() - pos.z(),
                         inter_p.x() - pos.x());  // between -M_PI ~ M_PI
        if (u < 0) u += 2 * M_PI;
        return Vector3f(u, t, 0);
    }
    double get_sphere_intersect(Ray ray, Vector3f o, double r) {
        Vector3f ro = o - ray.getOrigin();
        double b = Vector3f::dot(ray.getDirection(), ro);
        double d = sqr(b) - ro.squaredLength() + sqr(r);
        if (d < 0)
            return -1;
        else
            d = sqrt(d);
        double t = b - d > eps ? b - d : b + d > eps ? b + d : -1;
        if (t < 0) return -1;
        return t;
    }
    virtual bool intersect(const Ray& ray, Hit& hit, double tmin) {
        std::pair<double, Vector3f> temp = intersect(ray);
        double t = temp.first;
        if (t >= tmin && t < hit.getT()) {
            Vector3f n = norm(temp.second);
            if (Vector3f::dot(n, -ray.getDirection()) < 0) n = -n;

            hit.set(t, material, n);

            if (material->texture.valid()) {
                Vector3f texCoord = change_for_bezier(temp.second);
                hit.setTexCoord(
                    Vector2f(texCoord.x() / 2 / M_PI, texCoord.y()));
            }

            return 1;
        }
        return 0;
    }
    virtual std::pair<double, Vector3f> intersect(Ray ray) {
        double final_dis = INF;
        // check for |dy|<eps
        if (std::abs(ray.getDirection().y()) < 5e-4) {  // y 轴偏移很小
            double dis_to_axis =
                (Vector3f(pos.x(), ray.getOrigin().y(), pos.z()) -
                 ray.getOrigin())
                    .length();  // 源点到 y 轴的距离
            double hit = ray.pointAtParameter(dis_to_axis).y();  // 挪 一段距离

            if (hit < pos.y() + eps || hit > pos.y() + curve.height - eps)
                return std::make_pair(INF, Vector3f());

            // solve function pos.y+y(t)=ray.o.y to get x(t)
            double t = solve_t(hit - pos.y());  // 得到 y 对应参数
            if (t < 0 || t > 1) return std::make_pair(INF, Vector3f());
            Vector3f loc = curve.getpos(t);
            double ft = pos.y() + loc.y() - hit;
            if (std::abs(ft) > eps) return std::make_pair(INF, Vector3f());
            // assume sphere (pos.x, pos.y + loc.y, pos.z) - loc.x
            final_dis = get_sphere_intersect(
                ray, Vector3f(pos.x(), pos.y() + loc.y(), pos.z()),
                loc.x());  // 圆求交
            if (final_dis < 0) return std::make_pair(INF, Vector3f());
            Vector3f inter_p = ray.pointAtParameter(final_dis);
            // printf("y %f small!!!",std::abs((inter_p - Vector3f(pos.x,
            // inter_p.y, pos.z)).len2() - sqr(loc.x)));
            if (std::abs((inter_p - Vector3f(pos.x(), inter_p.y(), pos.z()))
                             .squaredLength() -
                         sqr(loc.x())) > 1e-1)
                return std::make_pair(INF, Vector3f());
            // second iteration, more accuracy
            // 再搞一次一样的
            hit = inter_p.y();
            if (hit < pos.y() + eps || hit > pos.y() + curve.height - eps)
                return std::make_pair(INF, Vector3f());
            t = solve_t(hit - pos.y());
            loc = curve.getpos(t);
            ft = pos.y() + loc.y() - hit;
            if (std::abs(ft) > eps) return std::make_pair(INF, Vector3f());
            final_dis = get_sphere_intersect(
                ray, Vector3f(pos.x(), pos.y() + loc.y(), pos.z()),
                loc.x());  // 求圆交点
            if (final_dis < 0) return std::make_pair(INF, Vector3f());
            inter_p = ray.pointAtParameter(final_dis);
            if (std::abs((inter_p - Vector3f(pos.x(), hit, pos.z()))
                             .squaredLength() -
                         sqr(loc.x())) > 1e-2)
                return std::make_pair(INF, Vector3f());
            // printf("---y %f small!!!",std::abs((inter_p - Vector3f(pos.x,
            // inter_p.y, pos.z)).len2() - sqr(loc.x)));
            return std::make_pair(final_dis, inter_p);
        }
        // printf("y big\n");
        // check for top circle: the plane is y=pos.y + curve.height
        // TODO
        // check for buttom circle: the plane is y=pos.y
        // TODO
        // normal case
        // calc ay^2+by+c
        double a = 0, b = 0, c = 0, t1, t2;
        // (xo-x'+xd/yd*(y-yo))^2 -> (t1+t2*y)^2 = x^2
        t1 = ray.getOrigin().x() - pos.x() -
             ray.getDirection().x() / ray.getDirection().y() *
                 ray.getOrigin().y();
        t2 = ray.getDirection().x() / ray.getDirection().y();
        a += t2 * t2;
        b += 2 * t1 * t2;
        c += t1 * t1;
        // (zo-z'+zd/yd*(y-yo))^2 -> (t1+t2*y)^2 = z^2
        t1 = ray.getOrigin().z() - pos.z() -
             ray.getDirection().z() / ray.getDirection().y() *
                 ray.getOrigin().y();
        t2 = ray.getDirection().z() / ray.getDirection().y();
        a += sqr(t2);
        b += 2 * t1 * t2;
        c += sqr(t1);
        // ay^2+by+c -> a'(y-b')^2+c' = x^2+z^2
        c = c - b * b / 4 / a;
        b = -b / 2 / a - pos.y();
        // printf("%lf %lf %lf\n",a,b,c);
        if (0 <= b && b <= curve.height && c > curve.max2 ||
            (b < 0 || b > curve.height) &&
                std::min(sqr(b), sqr(curve.height - b)) * a + c >
                    curve.max2)  // no intersect 无解
            return std::make_pair(INF, Vector3f());
        // double M_PIck[20] = {0, 0, 1}; int tot = 2;
        // for (double _ = 0; _ <= 1; _ += 0.1)
        // {
        // 	double t_M_PIck = newton2(_, a, b, c);
        // 	if (0 <= t_M_PIck && t_M_PIck <= 1)
        // 	{
        // 		bool flag = 1;
        // 		for (int j = 1; j <= tot; ++j)
        // 			if (std::abs(t_M_PIck - M_PIck[j]) < eps)
        // 				flag = 0;
        // 		if (flag)
        // 			M_PIck[++tot] = t_M_PIck;
        // 	}
        // }
        // std::sort(M_PIck + 1, M_PIck + 1 + tot);
        // for (int j = 1; j < tot; ++j)
        // 	if (getft(M_PIck[j], a, b, c) * getft(M_PIck[j + 1], a, b, c) <=
        // 0) 		check(M_PIck[j], M_PIck[j+1], (M_PIck[j] + M_PIck[j +
        // 1]) * .5, ray, a, b, c, final_dis);
        for (int ind = 0; ind <= curve.num; ++ind) {
            // y = curve.ckpt[ind] ~ curve.ckpt[ind+1]
            // calc min(a(y-b)^2+c)
            // double lower;
            // if (curve.data[ind].y0 <= b && b <= curve.data[ind].y1)
            // 	lower = c;
            // else
            // 	lower = a * std::min(sqr(curve.data[ind].y0 - b),
            // sqr(curve.data[ind].y1 - b)) + c;
            double t0 = curve.data[ind].t0, t1 = curve.data[ind].t1;
            // if (t0 > eps) t0 += erand48(mess) * .01;
            // if (t1 < 1 - eps) t1 -= erand48(mess) * .01;
            // if (lower <= curve.data[ind].width2)
            {
                check(t0, t1, (t0 + t1 + t0) / 3, ray, a, b, c, final_dis);
                check(t0, t1, (t1 + t0 + t1) / 3, ray, a, b, c, final_dis);
            }
        }
        if (final_dis < INF / 2)
            return std::make_pair(final_dis, ray.pointAtParameter(final_dis));
        else
            return std::make_pair(INF, Vector3f());
    }
    bool check(double low, double upp, double init, Ray ray, double a, double b,
               double c, double& final_dis) {
        double t = newton(init, a, b, c, low, upp);
        if (t <= 0 || t >= 1) return false;
        Vector3f loc = curve.getpos(t);
        double x = loc.x(), y = loc.y();
        double ft = x - sqrt(a * sqr(y - b) + c);
        if (std::abs(ft) > eps) return false;
        // calc t for ray
        double dis =
            (pos.y() + y - ray.getOrigin().y()) / ray.getDirection().y();
        if (dis < eps) return false;
        Vector3f inter_p = ray.pointAtParameter(dis);
        if (std::abs((Vector3f(pos.x(), pos.y() + y, pos.z()) - inter_p)
                         .squaredLength() -
                     x * x) > eps)
            return false;
        if (dis < final_dis) {
            final_dis = dis;
            // printf("%lf %lf %lf %lf\n",t,x , sqrt(a * sqr(y - b) + c), x -
            // sqrt(a * sqr(y - b) + c));
            return true;
        }
        return false;
    }
    double getft(double t, double a, double b, double c) {
        if (t < 0) t = eps;
        if (t > 1) t = 1 - eps;
        Vector3f loc = curve.getpos(t);
        double x = loc.x(), y = loc.y();
        return x - sqrt(a * sqr(y - b) + c);
    }
    double newton(double t, double a, double b, double c, double low = eps,
                  double upp = 1 - eps) {
        // solve sqrt(a(y(t)+pos.y-b)^2+c)=x(t)
        // f(t) = x(t) - sqrt(a(y(t)+pos.y-b)^2+c)
        // f'(t) = x'(t) - a(y(t)+pos.y-b)*y'(t) / sqrt(...)
        // if t is not in [0, 1] then assume f(t) is a linear function
        double ft, dft, x, y, dx, dy, sq;
        Vector3f loc, dir;
        for (int i = 10; i--;) {
            if (t < 0) t = low;
            if (t > 1) t = upp;
            loc = curve.getpos(t), dir = curve.getdir(t);
            x = loc.x(), dx = dir.x();
            y = loc.y(), dy = dir.y();
            // printf("%lf %lf %lf\n",t,x,y);
            sq = sqrt(a * sqr(y - b) + c);
            ft = x - sq;
            dft = dx - a * (y - b) * dy / sq;
            if (std::abs(ft) < eps) return t;
            t -= ft / dft;
        }
        return -1;
    }
    double newton2(double t, double a, double b, double c) {
        double dft, ddft, y, dx, dy, ddx, ddy, sq;
        Vector3f loc, dir, dir2;
        for (int i = 5; i--;) {
            if (t < 0) t = eps;
            if (t > 1) t = 1 - eps;
            loc = curve.getpos(t), dir = curve.getdir(t),
            dir2 = curve.getdir2(t);
            y = loc.y(), dx = dir.x(), dy = dir.y();
            ddx = dir2.x(), ddy = dir2.y();
            sq = sqrt(a * sqr(y - b) + c);
            dft = dx - a * (y - b) * dy / sq;
            ddft = ddx - a * ((y - b) * ddy + sqr(dy)) / sq +
                   sqr(a * (y - b) * dy) / sq / sq / sq;
            if (std::abs(dft) < eps) return t;
            t -= dft / ddft;
        }
        return -1;
    }
    virtual std::pair<Vector3f, Vector3f> aabb() {
        return std::make_pair(
            Vector3f(pos.x() - curve.max, pos.y(), pos.z() - curve.max),
            Vector3f(pos.x() + curve.max, pos.y() + curve.height,
                     pos.z() + curve.max));
    }
    virtual Vector3f norm(Vector3f p) {
        Vector3f tmp = change_for_bezier(p);
        Vector3f dir = curve.getdir(tmp.y());
        Vector3f d_surface =
            Vector3f(cos(tmp.x()), dir.y() / dir.x(), sin(tmp.x()));
        Vector3f d_circ = Vector3f(-sin(tmp.x()), 0, cos(tmp.x()));
        return Vector3f::cross(d_circ, d_surface).normalized();
    }
};

#endif  // BEZIER_H