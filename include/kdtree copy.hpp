#ifndef KDTREE_H
#define KDTREE_H

#include <aabb.hpp>
#include <algorithm>

Vector3f min(const Vector3f &x, const Vector3f &y) {
    return Vector3f(std::min(x.x(), y.x()), std::min(x.y(), y.y()),
                    std::min(x.z(), y.z()));
}
Vector3f max(const Vector3f &x, const Vector3f &y) {
    return Vector3f(std::max(x.x(), y.x()), std::max(x.y(), y.y()),
                    std::max(x.z(), y.z()));
}

class KDTree {
   private:
    struct Node {
        AABB box;
        Node *ls = nullptr, *rs = nullptr;
    };

    int dim;
    Node *root;

   public:
    KDTree(AABB *a[], const int &n) { root = build(a, n); }
    ~KDTree() { destroy(root); }

    bool comp(AABB *x, AABB *y) {  // 按照第 dim 维中点划分
        return x->pmid[dim] < y->pmid[dim];
    }
    Node *newNode() { return new Node; }
    Node *build(AABB *a[], const int &n) {
        if (n == 1) {
            Node *u = newNode();
            u->box = *a[0];
            return;
        }

        int mid = n / 2;
        dim = (dim + 1) % 3;
        std::nth_element(a, a + mid, a + n, comp);

        Node *u = newNode();
        u->ls = build(a, mid);
        u->rs = build(a + mid, n - mid);
        u->box = AABB(min(u->ls->box.pmin, u->rs->box.pmin),
                      max(u->ls->box.pmax, u->rs->box.pmax));
        return u;
    }
    bool intersect(Node *u, const Ray &r, Hit &h, const float &tmin) {
        if (!u->box.intersect(r)) return 0;
        if (!u->ls) return u->box.obj->intersect(r, h, tmin);

        return intersect(u->ls, r, h, tmin) | intersect(u->rs, r, h, tmin);
    }
    void destroy(Node *u) {
        if (!u->ls) {
            delete u;
            return;
        }

        destroy(u->ls);
        destroy(u->rs);
        delete u;
    }
};

#endif  // KDTREE_H