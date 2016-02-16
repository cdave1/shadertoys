// Via: https://www.shadertoy.com/view/XdB3Ww
//
//Find roots using Cardano's method. http://en.wikipedia.org/wiki/Cubic_function#Cardano.27s_method
#define EPSILON 1e-8
#define PI 3.14159265358979


// Test if point p crosses line (a, b), returns sign of result
float testCross(vec2 a, vec2 b, vec2 p) {
    return sign((b.y-a.y) * (p.x-a.x) - (b.x-a.x) * (p.y-a.y));
}

// Determine which side we're on (using barycentric parameterization)
float signBezier(vec2 A, vec2 B, vec2 C, vec2 p)
{
    vec2 a = C - A, b = B - A, c = p - A;
    vec2 bary = vec2(c.x*b.y-b.x*c.y,a.x*c.y-c.x*a.y) / (a.x*b.y-b.x*a.y);
    vec2 d = vec2(bary.y * 0.5, 0.0) + 1.0 - bary.x - bary.y;
    return mix(sign(d.x * d.x - d.y), mix(-1.0, 1.0,
        step(testCross(A, B, p) * testCross(B, C, p), 0.0)),
        step((d.x - d.y), 0.0)) * testCross(A, C, B);
}

// Solve cubic equation for roots
vec3 solveCubic(float a, float b, float c)
{
    float p = b - a*a / 3.0, p3 = p*p*p;
    float q = a * (2.0*a*a - 9.0*b) / 27.0 + c;
    float d = q*q + 4.0*p3 / 27.0;
    float offset = -a / 3.0;
    if(d >= 0.0) {
        float z = sqrt(d);
        vec2 x = (vec2(z, -z) - q) / 2.0;
        vec2 uv = sign(x)*pow(abs(x), vec2(1.0/3.0));
        return vec3(offset + uv.x + uv.y);
    }
    float v = acos(-sqrt(-27.0 / p3) * q / 2.0) / 3.0;
    float m = cos(v), n = sin(v)*1.732050808;
    return vec3(m + m, -n - m, n - m) * sqrt(-p / 3.0) + offset;
}

// Find the signed distance from a point to a bezier curve
float sdBezier(vec2 A, vec2 B, vec2 C, vec2 p)
{
    B = mix(B + vec2(1e-4), B, abs(sign(B * 2.0 - A - C)));
    vec2 a = B - A, b = A - B * 2.0 + C, c = a * 2.0, d = A - p;
    vec3 k = vec3(3.*dot(a,b),2.*dot(a,a)+dot(d,b),dot(d,a)) / dot(b,b);
    vec3 t = clamp(solveCubic(k.x, k.y, k.z), 0.0, 1.0);
    vec2 pos = A + (c + b*t.x)*t.x;
    float dis = length(pos - p);
    pos = A + (c + b*t.y)*t.y;
    dis = min(dis, length(pos - p));
    pos = A + (c + b*t.z)*t.z;
    dis = min(dis, length(pos - p));
    return dis * signBezier(A, B, C, p);
}


int findRoots(float a, float b, float c, float d, out float r[3])
{
	if (abs(a) > EPSILON)
	{
		float z = 1.0/a;
		float d3 = 1.0/3.0;
		float d27 = 1.0/27.0;
		a = b*z;
		b = c*z;
		c = d*z;
		float p = b-a*a*d3;
		float q = a*(2.0*a*a-9.0*b)*d27+c;
		float ppp = p*p*p;
		float D = q*q+4.0*ppp*d27;
		float delta = -a*d3;
		if (D > EPSILON)
		{
			z = sqrt(D);
			float u = (-q+z)*0.5;
			float v = (-q-z)*0.5;
			u = sign(u)*pow(abs(u),d3);
			v = sign(v)*pow(abs(v),d3);
			r[0] = u+v+delta;
			return 1;
		}
		else if (D < -EPSILON)
		{
			float u = sqrt(-p*d3)*2.0;
			float v = acos(-sqrt(-27.0/ppp)*q*0.5)*d3;
			r[0] = u*cos(v)+delta;
			r[1] = u*cos(v+2.0*PI*d3)+delta;
			r[2] = u*cos(v+4.0*PI*d3)+delta;
			return 3;
		}
		else
		{
			q = sign(q)*pow(abs(q)*0.5,d3);
			r[0] = 2.0*-q+delta;
			r[1] = q+delta;
			return 2;
		}
	}
	else
	{
		if (abs(b) <= EPSILON && abs(c) > EPSILON)
		{
			r[0] = -d/c;
			return 1;
		}
		else
		{
			float D = c*c-4.0*b*d;
			float z = 1.0/(2.0*b);
			if (D > EPSILON)
			{
				D = sqrt(D);
				r[0] = (-c-D)*z;
				r[1] = (-c+D)*z;
				return 2;
			}
			else if (D > -EPSILON)
			{
				r[0] = -c*z;
				return 1;
			}
		}
	}
	return 0;
}

vec2 getPositionOnBezierCurve(float t, vec2 p0, vec2 p1, vec2 p2)
{
	float fOneMinusT = 1.0-t;
	return fOneMinusT*fOneMinusT*p0+2.0*t*fOneMinusT*p1+t*t*p2;
}

// How to resolve the equation below can be seen on this image.
// http://www.perbloksgaard.dk/research/DistanceToQuadraticBezier.jpg
float calculateDistanceToQuadraticBezier(vec2 p, vec2 p0, vec2 p1, vec2 p2)
{
	vec2 A = p1-p0;
	vec2 B = p2-p1-A;
	vec2 C = p-p0;
	float a = -dot(B,B);
	float b = -3.0*dot(A,B);
	float c = dot(C,B)-2.0*dot(A,A);
	float d = dot(C,A);
	float r[3];
	findRoots(a,b,c,d,r);
	float dist = distance(p,getPositionOnBezierCurve(clamp(r[0],0.0,1.0),p0,p1,p2));
	dist = min(dist, distance(p,getPositionOnBezierCurve(clamp(r[1],0.0,1.0),p0,p1,p2)));
	dist = min(dist, distance(p,getPositionOnBezierCurve(clamp(r[2],0.0,1.0),p0,p1,p2)));
	return dist;
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float aspectRatio = iResolution.x / iResolution.y;
    vec2 percent = ((fragCoord.xy / iResolution.xy) - vec2(0.25,0.5));
    percent.x *= aspectRatio;

    vec2 mouse = (iMouse.xy / iResolution.xy);// - vec2(0.25,0.5);

    vec2 A = vec2(0.1,0.1);
    vec2 B = length(iMouse.xy) > 0.0 ? mouse : vec2(0.5, 0.9);
    vec2 C = vec2(0.9,0.5);
    vec2 D = vec2(0.9, 0.9);
    vec2 E = vec2(0.75, 0.75);

    float d1, d2, d3;

    if (fract(iGlobalTime) > 0.5) {
        vec2 p = fragCoord.xy / iResolution.xy; //(2.0*fragCoord.xy-iResolution.xy)/iResolution.y;
    	d1 = sdBezier(A, B, C, p);
    	d2 = sdBezier(C, B, D, p);
    	d3 = sdBezier(E, B, A, p);

        float d = 1.0 - max(d1, max(d2, d3));
    	d = min(d, step(d, 0.90));

        fragColor = vec4(d,d,d,1.0);// - sign(d)*vec4(0.1,0.4,0.7,1.0);
    	//fragColor *= (1.0 - exp(-4.0*abs(d))) * (0.8 + 0.2*cos(140.*d));
    	//fragColor = mix(fragColor, vec4(1.0), 1.0-smoothstep(0.0,0.02,abs(d)) );

    } else {

    	A.xy *= iResolution.xy;
    	B.xy *= iResolution.xy;
    	C.xy *= iResolution.xy;
    	D.xy *= iResolution.xy;
    	E.xy *= iResolution.xy;

    	d1 = calculateDistanceToQuadraticBezier(fragCoord.xy, A, B, C);
    	d2 = calculateDistanceToQuadraticBezier(fragCoord.xy, C, B, D);
   		d3 = calculateDistanceToQuadraticBezier(fragCoord.xy, E, B, A);

        d1 = 1.0 - (d1 * 0.01);
    	d2 = 1.0 - (d2 * 0.01);
    	d3 = 1.0 - (d3 * 0.01);

    	float d = 0.0;

    	d = max(d1, max(d2, d3));
    	d = min(d, step(d, 0.90));
    	fragColor = vec4(d,d,d,1.0);
    }
}
