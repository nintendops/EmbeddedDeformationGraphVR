attribute vec3 inputPosition;
attribute vec2 inputTexCoord;
attribute vec3 inputNormal;

uniform mat4 projection, modelview, normalMat;
uniform int mode;

varying vec4 forFragColor;

const vec3 lightPos = vec3(1.0, 1.0, 1.0);
const vec3 diffuseColor = vec3(0.5, 0.0, 0.0);
const vec3 specColor = vec3(1.0, 1.0, 1.0);

void main(){
  gl_Position = projection * modelview * vec4(inputPosition, 1.0);

  // all following gemetric computations are performed in the
  // camera coordinate system (aka eye coordinates)
  vec3 normal = vec3(normalMat * vec4(inputNormal, 0.0));
  vec4 vertPos4 = modelview * vec4(inputPosition, 1.0);
  vec3 vertPos = vec3(vertPos4) / vertPos4.w;
  vec3 lightDir = normalize(lightPos - vertPos);
  vec3 reflectDir = reflect(-lightDir, normal);
  vec3 viewDir = normalize(-vertPos);

  float lambertian = max(dot(lightDir,normal), 0.0);
  float specular = 0.0;
  
  if(lambertian > 0.0) {
    float specAngle = max(dot(reflectDir, viewDir), 0.0);
    specular = pow(specAngle, 4.0);

    // the exponent controls the shininess (try mode 2)
    if(mode == 2)  specular = pow(specAngle, 16.0);
       
    // according to the rendering equation we would need to multiply
    // with the the "lambertian", but this has little visual effect
    if(mode == 3) specular *= lambertian;
    // switch to mode 4 to turn off the specular component
    if(mode == 4) specular *= 0.0;
  }
  
  forFragColor = vec4(lambertian*diffuseColor + specular*specColor, 1.0);
}
          