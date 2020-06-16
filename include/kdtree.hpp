#ifndef KDTREE_H
#define KDTREE_H

#include <aabb.hpp>
#include <algorithm>
#include <vector>

int dim = 0;
bool comp(AABB x, AABB y) {  // 按照第 dim 维中点划分
    return x.pmid[dim] < y.pmid[dim];
}

class KDTree {
   private:
    struct Node {
        AABB box;
        Node *ls = nullptr, *rs = nullptr;
    };

    // static int dim;
    Node *root = nullptr;

   public:
    KDTree(std::vector<AABB>::iterator a, const int &n) {
        if (n) root = build(a, n);
    }
    ~KDTree() {
        if (root) destroy(root);
    }

    Vector3f min(const Vector3f &x, const Vector3f &y) {
        return Vector3f(std::min(x.x(), y.x()), std::min(x.y(), y.y()),
                        std::min(x.z(), y.z()));
    }
    Vector3f max(const Vector3f &x, const Vector3f &y) {
        return Vector3f(std::max(x.x(), y.x()), std::max(x.y(), y.y()),
                        std::max(x.z(), y.z()));
    }

    Node *newNode() { return new Node; }
    Node *build(std::vector<AABB>::iterator a, const int &n) {
        if (n == 1) {
            Node *u = newNode();
            u->box = *a;
            return u;
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
    bool intersect(const Ray &r, Hit &h, const double &tmin) const {
        if (root)
            return intersect(root, r, h, tmin);
        else
            return 0;
    }
    bool intersect(Node *u, const Ray &r, Hit &h, const double &tmin) const {
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