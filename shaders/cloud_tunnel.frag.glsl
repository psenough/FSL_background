#version 300 es
precision mediump float;
in vec2 vUV;
out vec4 fragColor;
uniform float iTime;

vec3 pal[6]= vec3[](vec3(1.0,1.0,1.0),
                        vec3(106.0/255.0, 169.0/255.0, 219.0/255.0),
                        vec3(30.0/255.0, 104.0/255.0, 175.0/255.0),
                        vec3(95.0/255.0, 152.0/255.0, 200.0/255.0),
                        vec3(57.0/255.0, 143.0/255.0, 200.0/255.0),
                        vec3(36.0/255.0, 52.0/255.0, 114.0/255.0));

vec2 iResolution = vec2(1.0,1.0);

vec3 palette_limiter (in vec3 albedo){
	float estimation_cutoff = 0.001;
	vec3 closest_color;
	float min_dist = 2.0;

	for (int i=0; i<6; i++ ){
		//float index = 1.000/(2.000*n)+float(i)/n;
		vec3 index_color = pal[i];
		float dist = length(index_color - albedo);
		if (dist < min_dist) {
			min_dist = dist;
			closest_color = index_color;
			if (min_dist < estimation_cutoff){
				return closest_color;
			}
		}
	}
	return closest_color;
}

void main() {
	//Raymarch depth
    float z,
    //Step distance
    d,
    //Raymarch iterator
    i,
    //Animation time
    t = iTime;

	vec4 O = vec4(0.0);

    //Clear fragColor and raymarch 100 steps
    for(O*=i; i++<1e2;
        //Coloring and brightness
        O += (1.+cos(i*.7+vec4(6,1,2,0)))/d/i)
    {
    
        //Sample point (from ray direction)
        vec3 p = z*normalize(vec3(vUV+vUV,0)-iResolution.xyx)+.1;
        
        //Polar coordinates
        p = vec3(atan(p.y,p.x)*2., p.z/3., length(p.xy)-6.);
        
        //Apply turbulence
        //https://mini.gmshaders.com/p/turbulence
        for(d=0.; d++<9.;)
            p += sin(p.yzx*d-t+.2*i)/d/3.0;
            
        //Distance to cylinder and waves
        z+=d=.2*length(vec4(p.z,.1*cos(p*3.)-.1));
    }
    //Tanh tonemap
    O=tanh(O*O/9e2);
    O=vec4(palette_limiter(O.rgb*1.5),O.a);


    fragColor = O;
}