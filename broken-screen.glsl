vec3 rand(vec2 coord) {
    vec3 a = fract(cos(coord.x*8.3e-3 + coord.y )*vec3(1.3e5, 4.7e5, 2.9e5));
    vec3 b = fract(sin(coord.x*8.3e-3 + coord.y )*vec3(8.1e5, 1.0e5, 0.1e5));
    vec3 c = mix(a, b, 0.5);
    return c;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord.xy / iResolution.xy;
    vec4 color = vec4(rand(uv *sin(10.0 * iGlobalTime)), 1.) * 2. - 1.;
    fragColor = color;
}
