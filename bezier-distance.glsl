// Via: https://www.shadertoy.com/view/XdB3Ww
//
//Find roots using Cardano's method. http://en.wikipedia.org/wiki/Cubic_function#Cardano.27s_method
#define EPSILON 1e-8
#define PI 3.14159265358979

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
    vec2 B = length(iMouse.xy) > 0.0 ? mouse : vec2(-0.3,0.2);
    vec2 C = vec2(0.9,0.5);
    vec2 D = vec2(0.9, 0.9);
    vec2 E = vec2(0.75, 0.75);

    A.xy *= iResolution.xy;
    B.xy *= iResolution.xy;
    C.xy *= iResolution.xy;
    D.xy *= iResolution.xy;
    E.xy *= iResolution.xy;

    float d1 = calculateDistanceToQuadraticBezier(fragCoord.xy, A, B, C);
    float d2 = calculateDistanceToQuadraticBezier(fragCoord.xy, C, B, D);
    float d3 = calculateDistanceToQuadraticBezier(fragCoord.xy, E, B, A);

    d1 = 1.0 - (d1 * 0.01);
    d2 = 1.0 - (d2 * 0.01);
    d3 = 1.0 - (d3 * 0.01);

    float d = 0.0;

    d = max(d1, max(d2, d3));

    if (d > 0.90)
        fragColor = vec4(0.0, 0.0, 0.0, 1.0);
    else
    	fragColor = vec4(d,d,d,1.0);
}
