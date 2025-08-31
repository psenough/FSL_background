#version 300 es
precision mediump float;
in vec2 vUV;
out vec4 fragColor;
uniform float iTime;


vec2 iResolution = vec2(1.0,1.0);

float sdCornerCircle(   in vec2 p ); // c(u) = sqrt(2-u^2)-1
float sdCornerParabola( in vec2 p ); // c(u) = (1-u^2)/2
float sdCornerCosine(   in vec2 p ); // c(u) = (2/PI)*cos(u*PI/2)
float sdCornerCubic(    in vec2 p ); // c(u) = (u^3-3u^2+2)/3

const float kT = 6.28318531;

float sdRoundBox( in vec2 p, in vec2 b, in vec4 r, in int type )
{
    // select corner radius
    r.xy = (p.x>0.0)?r.xy : r.zw;
    r.x  = (p.y>0.0)?r.x  : r.y;
    // box coordinates
    vec2 q = abs(p)-b+r.x;
    // distance to sides
    if( min(q.x,q.y)<0.0 ) return max(q.x,q.y)-r.x;
    // rotate 45 degrees, offset by r and scale by r*sqrt(0.5) to canonical corner coordinates
    vec2 uv = vec2( abs(q.x-q.y), q.x+q.y-r.x )/r.x;
    // compute distance to corner shape
    float d;
         if( type==0 ) d = sdCornerCircle( uv );
    else if( type==1 ) d = sdCornerParabola( uv );
    else if( type==2 ) d = sdCornerCosine( uv );
    else if( type==3 ) d = sdCornerCubic( uv );
    // undo scale
    return d * r.x*sqrt(0.5);
}

float sdCornerCircle( in vec2 p )
{
    return length(p-vec2(0.0,-1.0)) - sqrt(2.0);
}

float sdCornerParabola( in vec2 p )
{
    // https://www.shadertoy.com/view/ws3GD7
    float y = (0.5 + p.y)*(2.0/3.0);
    float h = p.x*p.x + y*y*y;
    float w = pow( p.x + sqrt(abs(h)), 1.0/3.0 ); // note I allow a tiny error in the very interior of the shape so that I don't have to branch into the 3 root solution
    float x = w - y/w;
    vec2  q = vec2(x,0.5*(1.0-x*x));
    return length(p-q)*sign(p.y-q.y);
}

float sdCornerCosine( in vec2 uv )
{
    // https://www.shadertoy.com/view/3t23WG
    uv *= (kT/4.0);

    float ta = 0.0, tb = kT/4.0;
    for( int i=0; i<8; i++ )
    {
        float t = 0.5*(ta+tb);
        float y = t-uv.x+(uv.y-cos(t))*sin(t);
        if( y<0.0 ) ta = t; else tb = t;
    }
    vec2  qa = vec2(ta,cos(ta)), qb = vec2(tb,cos(tb));
    vec2  pa = uv-qa, di = qb-qa;
    float h = clamp( dot(pa,di)/dot(di,di), 0.0, 1.0 );
    return length(pa-di*h) * sign(pa.y*di.x-pa.x*di.y) * (4.0/kT);
}

float sdCornerCubic( in vec2 uv )
{
    float ta = 0.0, tb = 1.0;
    for( int i=0; i<12; i++ )
    {
        float t = 0.5*(ta+tb);
        float c = (t*t*(t-3.0)+2.0)/3.0;
        float dc = t*(t-2.0);
        float y = (uv.x-t) + (uv.y-c)*dc;
        if( y>0.0 ) ta = t; else tb = t;
    }
    vec2  qa = vec2(ta,(ta*ta*(ta-3.0)+2.0)/3.0);
    vec2  qb = vec2(tb,(tb*tb*(tb-3.0)+2.0)/3.0);
    vec2  pa = uv-qa, di = qb-qa;
    float h = clamp( dot(pa,di)/dot(di,di),0.0,1.0 );
    return length(pa-di*h) * sign(pa.y*di.x-pa.x*di.y);
}

float approx_sdSuperEllipse( vec2 p, vec2 b, vec4 r )
{
    // select corner radius
    r.xy = (p.x>0.0)?r.xy : r.zw;
    r.x  = (p.y>0.0)?r.x  : r.y;
    
    float n = r.x;
    
    p = abs(p);
    
    #if 0
        // reall bad, cheap linearliation of the basic implicit formula
        n = 2.0/n;
        float w = pow(p.x/b.x,n) + pow(p.y/b.y,n);
        float kb = 2.0*n - 2.0;
        float ka = 1.0-1.0/n;
        float kc = 2.0*n;
        return (w-pow(w,ka)) * inversesqrt( pow(p.x,kb)/pow(b.x,kc) + pow(p.y,kb)/pow(b.y,kc) );
    #else
        // bisection root finder for distance minimizer
        float xa = 0.0, xb = kT/4.0;
        if( p.x-b.x>p.y-b.y ) xb*=0.5; // hack for interior distances, still wrong
        for( int i=0; i<12; i++ ) 
        {
            float x = 0.5*(xa+xb);
            float c = cos(x);
            float s = sin(x);
            float cn = pow(c,n);
            float sn = pow(s,n);
            float y = (p.x-b.x*cn)*b.x*cn*s*s - (p.y-b.y*sn)*b.y*sn*c*c;
            if( y<0.0 ) xa = x; else xb = x;
        }
        vec2  qa = b*pow(vec2(cos(xa),sin(xa)),vec2(n));
        vec2  qb = b*pow(vec2(cos(xb),sin(xb)),vec2(n));
        vec2  pa = p-qa, ba = qb-qa;
        float h = clamp( dot(pa,ba)/dot(ba,ba),0.0,1.0);
        return length( pa - ba*h ) * sign(pa.x*ba.y-pa.y*ba.x);
    #endif
}


float round_merge(float shape1, float shape2, float radius){
    vec2 intersectionSpace = vec2(shape1 - radius, shape2 - radius);
    intersectionSpace = min(intersectionSpace, 0.0);
    float insideDistance = -length(intersectionSpace);
    float simpleUnion = min(shape1, shape2);
    float outsideDistance = max(simpleUnion, radius);
    return  insideDistance + outsideDistance;
}

void main() {
	// coordinates
	vec2 p = (2.0*vUV-iResolution.xy)/iResolution.y*3.0;

    // darker cloud
    vec2 pd = p + vec2(-1.24+sin(iTime/8.0)*.6,2.12);
    float d = sdRoundBox( pd,                  vec2(1.2, 0.1), vec4(0.1,0.1,0.1,0.1), 1 );
    float d2 = sdRoundBox( pd+vec2(-1.0, 0.2), vec2(0.8,0.1), vec4(0.1,0.1,0.1,0.1), 1 );
    float d3 = sdRoundBox( pd+vec2(-.4, 0.4), vec2(1.0,0.1), vec4(0.1,0.1,0.1,0.1), 1 );
    float d4 = sdRoundBox( pd+vec2(-.2, -0.2), vec2(.6,0.2), vec4(0.2,0.2,0.2,0.2), 1 );
    float res = round_merge(d, d2, 0.1);
    res=round_merge(res,d3,0.1);
    res=round_merge(res,d4,0.1);
    
    
    // lighter cloud
    vec2 pc = p + vec2(2.28+sin(iTime/6.0+100.0)*.4,-2.42);
    float dc = sdRoundBox( pc+vec2(-.4,-0.1),    vec2(.6, 0.1), vec4(0.1,0.1,0.1,0.1), 1 );
    float dc2 = sdRoundBox( pc+vec2(-1.8, 0.), vec2(0.8,0.1), vec4(0.1,0.1,0.1,0.1), 1 );
    float dc3 = sdRoundBox( pc+vec2(-1.2, 0.4), vec2(1.2,0.1), vec4(0.1,0.1,0.1,0.1), 1 );
    float dc4 = sdRoundBox( pc+vec2(-1.3, 0.2), vec2(.8,0.2), vec4(0.2,0.2,0.1,0.2), 1 );
    float dc5 = sdRoundBox( pc+vec2(-.6, 0.6), vec2(1.,0.2), vec4(0.2,0.2,0.2,0.2), 1 );
    
    float resc = round_merge(dc, dc2, 0.1);
    resc=round_merge(resc,dc3,0.1);
    resc=round_merge(resc,dc4,0.1);
    resc=round_merge(resc,dc5,0.1);

    vec3 colc = (resc>0.0) ? vec3(34./255.,62./255.,120./255.) : vec3(96./255.,165./255.,215./255.);
    colc = mix(colc,vec3(34./255.,62./255.,120./255.),1.0-smoothstep(0.0,0.01,abs(resc)));

    vec3 col = (res>0.0) ? colc : vec3(57./255.,143./255.,197./255.);
    col = mix(col,colc,1.0-smoothstep(0.0,0.01,abs(res)));
	
    //col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,0.01,abs(res)) );

    // output
	fragColor = vec4(col,1.0);
    
}