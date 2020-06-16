// Normal Distribution Function
float NDF(const vec3 & h, float roughness) {
    //  GGX / Trowbridge-Reitz
    
    float alpha = roughness * roughness;
    float alpha2 = alpha * alpha;
    float NoH = h.z;
    return alpha2 / (PI*pow(NoH*NoH*(alpha2 - 1) + 1, 2));
}

// Fresnel
vec3 Fr(const vec3 & wi, const vec3 & h, const vec3 & albedo, float metallic) {
    // Schlick’s approximation
    // use a Spherical Gaussian approximation to replace the power.
    // slightly more efficient to calculate and the difference is imperceptible
    
    vec3 F0 = mix(vec3(0.04f), albedo, metallic);
    float HoWi = dot(h, wi);
    return F0 + (vec3(1.0f) - F0) * pow(2.0f, (-5.55473f * HoWi - 6.98316f) * HoWi);
}

// Geometry Function
float G(const vec3 & wo, const vec3 & wi, float roughness) {
    // Schlick, remap roughness and k

    // k = alpha / 2
    // direct light: alpha = pow( (roughness + 1) / 2, 2)
    // IBL(image base lighting) : alpha = pow( roughness, 2)

    if (wo.z <= 0 || wi.z <= 0)
        return 0;

    float k = pow(roughness + 1, 2) / 8.f;
    float G1_wo = wo.z / (wo.z*(1 - k) + k);
    float G1_wi = wi.z / (wi.z*(1 - k) + k);
    return G1_wo * G1_wi;
}

// cook torrance BRDF
vec3 CT_BRDF(const vec3 & wo, const vec3 & wi, const vec3 & albedo,
    float metallic, float roughness)
{
    vec3 h = normalize(wo + wi);
    return NDF(h, roughness)*Fr(wi, h, albedo, metallic)*G(wo, wi, roughness) / (4 * wo.z*wi.z);
}

// 采样 BRDF
// pd 是 probability density
vec3 Sample_f(const vec3 & wo,
    const vec3 & albedo, float roughness, float metallic, float ambientOcclusion,
    vec3 & wi, float & pd)
{
    // 根据 NDF 采样 h
    float Xi1 = Rand_F();
    float Xi2 = Rand_F();
    float alpha = roughness * roughness;
    float cosTheta2 = (1 - Xi1) / (Xi1*(alpha*alpha - 1) + 1);
    float cosTheta = sqrt(cosTheta2);
    float sinTheta = sqrt(1 - cosTheta2);
    float phi = 2 * PI * Xi2;
    vec3 h(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta);
    wi = reflect(-wo, h);
    if (wi.z <= 0) {
        pd = 0;
        return vec3(0);
    }
    pd = NDF(h, roughness) / 4.0f;

    vec3 diffuse = albedo / PI;
    return ambientOcclusion * ((1 - metallic)*diffuse +
        CT_BRDF(wo, wi, albedo, metallic, roughness));
}

// 概率密度函数 probability density function
float PDF(const vec3 & wo, const vec3 & wi) {
    vec3 h = normalize(wo + wi);
    return NDF(h, roughness) / 4.0f;
}

// BRDF
vec3 F(const vec3 & wo, const vec3 & wi, const vec3 & albedo,
    float roughness, float metallic, float ambientOcclusion)
{
    auto diffuse = albedo / PI;
    return ambientOcclusion * ((1 - metallic)*diffuse +
        CT_BRDF(wo, wi, albedo, metallic, roughness));
}
