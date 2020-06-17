class BezierCurve2D
{// f(y)=x, y goes up?
public:
	ld *dx, *dy, max, height, max2, r, num;
	int n;
	struct D{
		ld t0, t1, width, y0, y1, width2;
	}data[20];
	// x(t) = \sum_{i=0}^n dx_i * t^i
	// y(t) = \sum_{i=0}^n dy_i * t^i
	BezierCurve2D(ld* px, ld* py, int n_, int num_, ld r_): num(num_), n(n_), r(r_) {
		dx = new ld[n];
		dy = new ld[n];
		assert(std::abs(py[0]) <= 1e-6);
		--n;
		// preproces
		for(int i = 0; i <= n; ++i)
		{
			dx[i] = px[0];
			dy[i] = py[0];
			for (int j = 0; j <= n - i; ++j)
			{
				px[j] = px[j + 1] - px[j];
				py[j] = py[j + 1] - py[j];
			}
		}
		ld n_down = 1, fac = 1, nxt = n;
		for (int i = 0; i <= n; ++i, --nxt)
		{
			fac = fac * (i == 0 ? 1 : i);
			dx[i] = dx[i] * n_down / fac;
			dy[i] = dy[i] * n_down / fac;
			n_down *= nxt;
		}
		max = 0;
		ld interval = 1. / (num - 1), c = 0;
		for (int cnt = 0; cnt <= num; c += interval, ++cnt)
		{
			data[cnt].width = 0;
			data[cnt].t0 = std::max(0., c - r);
			data[cnt].t1 = std::min(1., c + r);
			data[cnt].y0 = getpos(data[cnt].t0).y;
			data[cnt].y1 = getpos(data[cnt].t1).y;
			for (ld t = data[cnt].t0; t <= data[cnt].t1; t += 0.00001)
			{
				P3 pos = getpos(t);
				if (data[cnt].width < pos.x)
					data[cnt].width = pos.x;
			}
			if (max < data[cnt].width)
				max = data[cnt].width;
			data[cnt].width += eps;
			data[cnt].width2 = sqr(data[cnt].width);
		}
		max += eps;
		max2 = max * max;
		height = getpos(1).y;
	}
	P3 getpos(ld t)
	{
		ld ans_x = 0, ans_y = 0, t_pow = 1;
		for (int i = 0; i <= n; ++i)
		{
			ans_x += dx[i] * t_pow;
			ans_y += dy[i] * t_pow;
			t_pow *= t;
		}
		return P3(ans_x, ans_y);
	}
	P3 getdir(ld t)
	{
		ld ans_x = 0, ans_y = 0, t_pow = 1;
		for(int i = 1; i <= n; ++i)
		{
			ans_x += dx[i] * i * t_pow;
			ans_y += dy[i] * i * t_pow;
			t_pow *= t;
		}
		return P3(ans_x, ans_y);
	}
	P3 getdir2(ld t)
	{
		ld ans_x = 0, ans_y = 0, t_pow = 1;
		for(int i = 2; i <= n; ++i)
		{
			ans_x += dx[i] * i * (i - 1) * t_pow;
			ans_y += dy[i] * i * (i - 1) * t_pow;
			t_pow *= t;
		}
		return P3(ans_x, ans_y);
	}
};


class BezierObject: public Object {
// the curve will rotate line (x=pos.x and z=pos.z) as pivot
public:
	BezierCurve2D curve;
	P3 pos; // the buttom center point
	BezierObject(P3 pos_, BezierCurve2D c_, Texture t):
		pos(pos_), curve(c_), Object(t) {}
	BezierObject(P3 pos_, BezierCurve2D c_, Refl_t refl, ld brdf = 1.5, P3 color = P3(), P3 emission = P3(), std::string tname = ""):
		pos(pos_), curve(c_), Object(refl, color, emission, brdf, tname) {}
	ld solve_t(ld yc) { // solve y(t)=yc
		// assert(0 <= yc && yc <= curve.height);
		ld t = .5, ft, dft;
		for (int i = 10; i--; )
		{
			if (t < 0) t = 0;
			else if (t > 1) t = 1;
			ft = curve.getpos(t).y - yc, dft = curve.getdir(t).y;
			if (std::abs(ft) < eps)
				return t;
			t -= ft / dft;
		}
		return -1;
	}
	virtual P3 change_for_bezier(P3 inter_p) {
		ld t = solve_t(inter_p.y - pos.y);
		ld u = atan2(inter_p.z - pos.z, inter_p.x - pos.x); // between -pi ~ pi
		if (u < 0)
			u += 2 * PI;
		return P3(u, t);
	}
	ld get_sphere_intersect(Ray ray, P3 o, ld r) {
		P3 ro = o - ray.o;
		ld b = ray.d.dot(ro);
		ld d = sqr(b) - ro.dot(ro) + sqr(r);
		if (d < 0) return -1;
		else d = sqrt(d);
		ld t = b - d > eps ? b - d : b + d > eps? b + d : -1;
		if (t < 0)
			return -1;
		return t;
	}
	virtual std::pair<ld, P3> intersect(Ray ray) {
		ld final_dis = INF;
		// check for |dy|<eps
		if (std::abs(ray.d.y) < 5e-4)
		{
			ld dis_to_axis = (P3(pos.x, ray.o.y, pos.z) - ray.o).len();
			ld hit = ray.get(dis_to_axis).y;
			if (hit < pos.y + eps || hit > pos.y + curve.height - eps)
				return std::make_pair(INF, P3());
			// solve function pos.y+y(t)=ray.o.y to get x(t)
			ld t = solve_t(hit - pos.y);
			if (t < 0 || t > 1)
				return std::make_pair(INF, P3());
			P3 loc = curve.getpos(t);
			ld ft = pos.y + loc.y - hit;
			if (std::abs(ft) > eps)
				return std::make_pair(INF, P3());
			// assume sphere (pos.x, pos.y + loc.y, pos.z) - loc.x
			final_dis = get_sphere_intersect(ray, P3(pos.x, pos.y + loc.y, pos.z), loc.x);
			if (final_dis < 0)
				return std::make_pair(INF, P3());
			P3 inter_p = ray.get(final_dis);
			// printf("y %f small!!!",std::abs((inter_p - P3(pos.x, inter_p.y, pos.z)).len2() - sqr(loc.x)));
			if (std::abs((inter_p - P3(pos.x, inter_p.y, pos.z)).len2() - sqr(loc.x)) > 1e-1)
				return std::make_pair(INF, P3());
			// second iteration, more accuracy
			hit = inter_p.y;
			if (hit < pos.y + eps || hit > pos.y + curve.height - eps)
				return std::make_pair(INF, P3());
			t = solve_t(hit - pos.y);
			loc = curve.getpos(t);
			ft = pos.y + loc.y - hit;
			if (std::abs(ft) > eps)
				return std::make_pair(INF, P3());
			final_dis = get_sphere_intersect(ray, P3(pos.x, pos.y + loc.y, pos.z), loc.x);
			if (final_dis < 0)
				return std::make_pair(INF, P3());
			inter_p = ray.get(final_dis);
			if (std::abs((inter_p - P3(pos.x, hit, pos.z)).len2() - sqr(loc.x)) > 1e-2)
				return std::make_pair(INF, P3());
			// printf("---y %f small!!!",std::abs((inter_p - P3(pos.x, inter_p.y, pos.z)).len2() - sqr(loc.x)));
			return std::make_pair(final_dis, inter_p);
		}
		// printf("y big\n");
		// check for top circle: the plane is y=pos.y + curve.height
		// TODO
		// check for buttom circle: the plane is y=pos.y
		// TODO
		// normal case
		// calc ay^2+by+c
		ld a = 0, b = 0, c = 0, t1, t2;
		// (xo-x'+xd/yd*(y-yo))^2 -> (t1+t2*y)^2
		t1 = ray.o.x - pos.x - ray.d.x / ray.d.y * ray.o.y;
		t2 = ray.d.x / ray.d.y;
		a += t2 * t2;
		b += 2 * t1 * t2;
		c += t1 * t1;
		// (zo-z'+zd/yd*(y-yo))^2 -> (t1+t2*y)^2
		t1 = ray.o.z - pos.z - ray.d.z / ray.d.y * ray.o.y;
		t2 = ray.d.z / ray.d.y;
		a += sqr(t2);
		b += 2 * t1 * t2;
		c += sqr(t1);
		// ay^2+by+c -> a'(y-b')^2+c'
		c = c - b * b / 4 / a;
		b = -b / 2 / a - pos.y;
		// printf("%lf %lf %lf\n",a,b,c);
		if (0 <= b && b <= curve.height && c > curve.max2
		 || (b < 0 || b > curve.height) && std::min(sqr(b), sqr(curve.height - b)) * a + c > curve.max2) // no intersect
			return std::make_pair(INF, P3());
		// ld pick[20] = {0, 0, 1}; int tot = 2;
		// for (ld _ = 0; _ <= 1; _ += 0.1)
		// {
		// 	ld t_pick = newton2(_, a, b, c);
		// 	if (0 <= t_pick && t_pick <= 1)
		// 	{
		// 		bool flag = 1;
		// 		for (int j = 1; j <= tot; ++j)
		// 			if (std::abs(t_pick - pick[j]) < eps)
		// 				flag = 0;
		// 		if (flag)
		// 			pick[++tot] = t_pick;
		// 	}
		// }
		// std::sort(pick + 1, pick + 1 + tot);
		// for (int j = 1; j < tot; ++j)
		// 	if (getft(pick[j], a, b, c) * getft(pick[j + 1], a, b, c) <= 0)
		// 		check(pick[j], pick[j+1], (pick[j] + pick[j + 1]) * .5, ray, a, b, c, final_dis);
		for(int ind = 0; ind <= curve.num; ++ind)
		{
			// y = curve.ckpt[ind] ~ curve.ckpt[ind+1]
			// calc min(a(y-b)^2+c)
			// ld lower;
			// if (curve.data[ind].y0 <= b && b <= curve.data[ind].y1)
			// 	lower = c;
			// else
			// 	lower = a * std::min(sqr(curve.data[ind].y0 - b), sqr(curve.data[ind].y1 - b)) + c;
			ld t0 = curve.data[ind].t0, t1 = curve.data[ind].t1;
			// if (t0 > eps) t0 += erand48(mess) * .01;
			// if (t1 < 1 - eps) t1 -= erand48(mess) * .01;
			// if (lower <= curve.data[ind].width2)
			{
				check(t0, t1, (t0 + t1 + t0) / 3, ray, a, b, c, final_dis);
				check(t0, t1, (t1 + t0 + t1) / 3, ray, a, b, c, final_dis);
			}
		}
		if (final_dis < INF / 2)
			return std::make_pair(final_dis, ray.get(final_dis));
		else
			return std::make_pair(INF, P3());
	}
	bool check(ld low, ld upp, ld init, Ray ray, ld a, ld b, ld c, ld&final_dis)
	{
		ld t = newton(init, a, b, c, low, upp);
		if (t <= 0 || t >= 1)
			return false;
		P3 loc = curve.getpos(t);
		ld x = loc.x, y = loc.y;
		ld ft = x - sqrt(a * sqr(y - b) + c);
		if (std::abs(ft) > eps)
			return false;
		// calc t for ray
		ld dis = (pos.y + y - ray.o.y) / ray.d.y;
		if (dis < eps)
			return false;
		P3 inter_p = ray.get(dis);
		if (std::abs((P3(pos.x, pos.y + y, pos.z) - inter_p).len2() - x * x) > eps)
			return false;
		if (dis < final_dis)
		{
			final_dis = dis;
			// printf("%lf %lf %lf %lf\n",t,x , sqrt(a * sqr(y - b) + c), x - sqrt(a * sqr(y - b) + c));
			return true;
		}
		return false;
	}
	ld getft(ld t, ld a, ld b, ld c)
	{
		if (t < 0) t = eps;
		if (t > 1) t = 1 - eps;
		P3 loc = curve.getpos(t);
		ld x = loc.x, y = loc.y;
		return x - sqrt(a * sqr(y - b) + c);
	}
	ld newton(ld t, ld a, ld b, ld c, ld low=eps, ld upp=1-eps)
	{
		// solve sqrt(a(y(t)+pos.y-b)^2+c)=x(t)
		// f(t) = x(t) - sqrt(a(y(t)+pos.y-b)^2+c)
		// f'(t) = x'(t) - a(y(t)+pos.y-b)*y'(t) / sqrt(...)
		// if t is not in [0, 1] then assume f(t) is a linear function
		ld ft, dft, x, y, dx, dy, sq;
		P3 loc, dir;
		for (int i = 10; i--; )
		{
			if (t < 0) t = low;
			if (t > 1) t = upp;
			loc = curve.getpos(t), dir = curve.getdir(t);
			x = loc.x, dx = dir.x;
			y = loc.y, dy = dir.y;
			// printf("%lf %lf %lf\n",t,x,y);
			sq = sqrt(a * sqr(y - b) + c);
			ft = x - sq;
			dft = dx - a * (y - b) * dy / sq;
			if (std::abs(ft) < eps)
				return t;
			t -= ft / dft;
		}
		return -1;
	}
	ld newton2(ld t, ld a, ld b, ld c)
	{
		ld dft, ddft, y, dx, dy, ddx, ddy, sq;
		P3 loc, dir, dir2;
		for (int i = 5; i--; )
		{
			if (t < 0) t = eps;
			if (t > 1) t = 1 - eps;
			loc = curve.getpos(t), dir = curve.getdir(t), dir2 = curve.getdir2(t);
			y = loc.y, dx = dir.x, dy = dir.y;
			ddx = dir2.x, ddy = dir2.y;
			sq = sqrt(a * sqr(y - b) + c);
			dft = dx - a * (y - b) * dy / sq;
			ddft = ddx - a * ((y - b) * ddy + sqr(dy)) / sq + sqr(a * (y - b) * dy) / sq / sq / sq;
			if (std::abs(dft) < eps)
				return t;
			t -= dft / ddft;
		}
		return -1;
	}
	virtual std::pair<P3, P3> aabb() {
		return std::make_pair(P3(pos.x - curve.max, pos.y, pos.z - curve.max), P3(pos.x + curve.max, pos.y + curve.height, pos.z + curve.max));
	}
	virtual P3 norm(P3 p) {
		P3 tmp = change_for_bezier(p);
		P3 dir = curve.getdir(tmp.y);
		P3 d_surface = P3(cos(tmp.x), dir.y / dir.x, sin(tmp.x));
		P3 d_circ = P3(-sin(tmp.x), 0, cos(tmp.x));
		return d_circ.cross(d_surface).norm();
	}
};