#ifndef VEC3_STRUCT
#define VEC3_STRUCT

#include <math.h>
#include <string>

#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

struct vec3 {
    float x, y, z;

    CUDA_HOSTDEV vec3 () : x(0.0f), y(0.0f), z(0.0f) {}

    CUDA_HOSTDEV vec3 (const vec3& v) : x(v.x), y(v.y), z(v.z) {}

    CUDA_HOSTDEV vec3 (float a, float b, float c) : x(a), y(b), z(c) {}

    CUDA_HOSTDEV void set(float a, float b, float c) {
        this->x = a;
        this->y = b;
        this->z = c;
    }
    std::string toString() const {
        return std::to_string(this->x) + "," + std::to_string(this->y) + "," + std::to_string(this->z);
    }
    CUDA_HOSTDEV void normalize() {
        float n = this->norm();
        this->x /= n;
        this->y /= n;
        this->z /= n;
    }
    CUDA_HOSTDEV inline vec3 operator+(const vec3& v) const {
        return vec3(this->x + v.x, this->y + v.y, this->z + v.z);
    }
    CUDA_HOSTDEV inline vec3 operator-(const vec3& v) const {
        return vec3(this->x - v.x, this->y - v.y, this->z - v.z);
    }
    CUDA_HOSTDEV inline vec3 operator*(const float& f) const {
        return vec3(this->x * f, this->y * f, this->z * f);
    }
    CUDA_HOSTDEV inline vec3 operator/(const float& f) const {
        return vec3(this->x / f, this->y / f, this->z / f);
    }
    CUDA_HOSTDEV inline vec3 operator+=(const float& a) {
        this->x += a;
        this->y += a;
        this->z += a;
        return *this;
    }
    CUDA_HOSTDEV inline vec3 operator-=(const float& a) {
        this->x -= a;
        this->y -= a;
        this->z -= a;
        return *this;
    }
    CUDA_HOSTDEV inline vec3 operator*=(const float& a) {
        this->x *= a;
        this->y *= a;
        this->z *= a;
        return *this;
    }
    CUDA_HOSTDEV inline vec3 operator/=(const float& a) {
        this->x /= a;
        this->y /= a;
        this->z /= a;
        return *this;
    }
    CUDA_HOSTDEV inline vec3 operator+=(const vec3& v) {
        this->x += v.x;
        this->y += v.y;
        this->z += v.z;
        return *this;
    }
    CUDA_HOSTDEV inline vec3 operator-=(const vec3& v) {
        this->x -= v.x;
        this->y -= v.y;
        this->z -= v.z;
        return *this;
    }
    CUDA_HOSTDEV inline vec3 operator*=(const vec3& v) {
        this->x *= v.x;
        this->y *= v.y;
        this->z *= v.z;
        return *this;
    }
    CUDA_HOSTDEV inline vec3 operator/=(const vec3& v) {
        this->x /= v.x;
        this->y /= v.y;
        this->z /= v.z;
        return *this;
    }
    CUDA_HOSTDEV inline float dot(const vec3& v) const {
        return this->x * v.x + this->y * v.y + this->z * v.z;
    }
    CUDA_HOSTDEV inline vec3 cross(const vec3& v) const {
        return vec3(this->y * v.z - this->z * v.y, this->z * v.x - this->x * v.z, this->x * v.y - this->y * v.x);
    }
    CUDA_HOSTDEV inline float norm() const {
        return sqrtf(this->x * this->x + this->y * this->y + this->z * this->z);
    }
    CUDA_HOSTDEV inline vec3 rotate(const float angle, const vec3& axis) const {
        // rotate vector about axis by angle
        float s = sinf(angle);
        float c = cosf(angle);
        // float d = this->dot(axis);
        // vec3 cr = axis.cross(*this);
        // return *this * c + cr * s + axis * d * (1.0f - c);
        return vec3(
            this->x * (c + (1.0f - c) * axis.x * axis.x) + this->y * ((1.0f - c) * axis.x * axis.y - s * axis.z) + this->z * ((1.0f - c) * axis.x * axis.z + s * axis.y),
            this->x * ((1.0f - c) * axis.y * axis.x + s * axis.z) + this->y * (c + (1.0f - c) * axis.y * axis.y) + this->z * ((1.0f - c) * axis.y * axis.z - s * axis.x),
            this->x * ((1.0f - c) * axis.z * axis.x - s * axis.y) + this->y * ((1.0f - c) * axis.z * axis.y + s * axis.x) + this->z * (c + (1.0f - c) * axis.z * axis.z)
        );
    }

};

#endif