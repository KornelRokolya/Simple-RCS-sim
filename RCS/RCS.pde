double t = 0;
double dt = 1e-2;
float heading = PI*1.1;
float pitch = -0.1;
float distance = 10;
float a = 1, b = 1, c = 1;
double R = 50;
double stdGravParam = 200;
PImage skybox, earth;
PFont font;
RigidBody body;
int viewMode = 0;
int iterationsPerGUI = 1;

void setup() {
  size(1500, 900, P3D);
  hint(ENABLE_DEPTH_MASK);
  skybox = loadImage("Skybox2.png");
  earth = loadImage("Earth.jpg");
  font = createFont("Arial", 20);
  
  body = new RigidBody();
  body.mass = 1;
  body.basis = Identity3x3();
  body.Theta0 = new Matrix3x3(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1);
  body.position = new Vector3(-80, 1, 1);
  body.velocity = getOrbitalVelocity(body.position);
  body.angularMoment = new Vector3();
  body.torque = new Vector3[0];
  body.torqueActive = new boolean[0];
  
  engines = new RCS_engine[] {
    new RCS_engine(-1,  0,  0,  a/2,    0,  c/2),
    new RCS_engine(-1,  0,  0,  a/2,  b/2,    0),
    new RCS_engine(-1,  0,  0,  a/2,    0, -c/2),
    new RCS_engine(-1,  0,  0,  a/2, -b/2,    0),
    new RCS_engine(-1,  0,  0,  a/2,    0,    0),
    new RCS_engine( 1,  0,  0, -a/2,    0,  c/2),
    new RCS_engine( 1,  0,  0, -a/2,  b/2,    0),
    new RCS_engine( 1,  0,  0, -a/2,    0, -c/2),
    new RCS_engine( 1,  0,  0, -a/2, -b/2,    0),
    new RCS_engine( 1,  0,  0, -a/2,    0,    0),
    new RCS_engine( 0, -1,  0,    0,  b/2,  c/2),
    new RCS_engine( 0, -1,  0,    0,  b/2, -c/2),
    new RCS_engine( 0, -1,  0,    0,  b/2,    0),
    new RCS_engine( 0,  1,  0,    0, -b/2,  c/2),
    new RCS_engine( 0,  1,  0,    0, -b/2, -c/2),
    new RCS_engine( 0,  1,  0,    0, -b/2,    0),
    new RCS_engine( 0,  0, -1,    0,    0,  c/2),
    new RCS_engine( 0,  0,  1,    0,    0, -c/2)
  };
}

void draw() {
  background(0);
  
  manageInput();
  for(int i = 0; i < iterationsPerGUI; i++) physicsUpdate();
  
  
  perspective(PI/3, 1.0*width/height, 0.01, 5000);
  lights();
  
  Matrix3x3 basis = body.basis;
  Vector3 position = body.position;
  Vector3 eyePosition = null;
  Vector3 center = null;
  Vector3 up = null;
  
  if(viewMode == 2) {
    eyePosition = lin(position, 1, basis.column0(), 0.51);
    center = lin(position, 1, basis.column0(), 1.51);
    up = basis.column2();
  }
  else if(viewMode == 1){
    eyePosition = new Vector3(
      cos(heading) * cos(pitch),
      sin(heading) * cos(pitch),
      sin(pitch)).mult(distance);
    center = new Vector3();
    up = new Vector3(
      -cos(heading)*sin(pitch),
      -sin(heading)*sin(pitch),
      cos(pitch));
  }
  else {
    eyePosition = position.add(new Vector3(
      cos(heading) * cos(pitch),
      sin(heading) * cos(pitch),
      sin(pitch)).mult(distance));
    center = position.Copy();
    up = new Vector3(
      -cos(heading)*sin(pitch),
      -sin(heading)*sin(pitch),
      cos(pitch));
  }
  
  camera(
      (float)eyePosition.x, (float)eyePosition.y, (float)eyePosition.z,
      (float)center.x, (float)center.y, (float)center.z,
      (float)up.x, (float)up.y, (float)up.z
    );
  
  hint(DISABLE_DEPTH_MASK);
  pushMatrix();
  translate((float)eyePosition.x, (float)eyePosition.y, (float)eyePosition.z);
  noStroke();
  fill(255, 255, 255);
  textureMode(IMAGE);
  textureWrap(REPEAT);
  PShape sky = createShape(BOX, 200);
  sky.setTexture(skybox);
  shape(sky);
  popMatrix();
  hint(ENABLE_DEPTH_MASK);
  
  
  // X tengely
  stroke(255,0,0);
  strokeWeight(2);
  beginShape(LINES);
  vertex(0,0,0);
  vertex(3,0,0);
  endShape();
  
  // Y tengely
  stroke(0,255,0);
  strokeWeight(2);
  beginShape(LINES);
  vertex(0,0,0);
  vertex(0,3,0);
  endShape();
  
  // Z tengely
  stroke(0,0,255);
  strokeWeight(2);
  beginShape(LINES);
  vertex(0,0,0);
  vertex(0,0,3);
  endShape();
  
  // Planet
  noStroke();
  PShape planet = createShape(SPHERE, (float)R);
  planet.setTexture(earth);
  pushMatrix();
  rotateX(PI/2);
  shape(planet);
  popMatrix();
  
  applyMatrix(
    (float)basis.a00, (float)basis.a01, (float)basis.a02, (float)body.position.x,
    (float)basis.a10, (float)basis.a11, (float)basis.a12, (float)body.position.y,
    (float)basis.a20, (float)basis.a21, (float)basis.a22, (float)body.position.z,
    0, 0, 0, 1
  );
  
  
  //shape(valami);
  stroke(255);
  fill(255,0,0);
  beginShape(QUAD);
  vertex(a/2, -b/2, -c/2);
  vertex(a/2, b/2, -c/2);
  vertex(a/2, b/2, c/2);
  vertex(a/2, -b/2, c/2);
  endShape(CLOSE);
  
  fill(128,0,0);
  beginShape(QUAD);
  vertex(-a/2, -b/2, -c/2);
  vertex(-a/2, b/2, -c/2);
  vertex(-a/2, b/2, c/2);
  vertex(-a/2, -b/2, c/2);
  endShape(CLOSE);
  
  fill(0,255,0);
  beginShape(QUAD);
  vertex(-a/2, b/2, -c/2);
  vertex(a/2, b/2, -c/2);
  vertex(a/2, b/2, c/2);
  vertex(-a/2, b/2, c/2);
  endShape(CLOSE);
  
  fill(0,128,0);
  beginShape(QUAD);
  vertex(-a/2, -b/2, -c/2);
  vertex(a/2, -b/2, -c/2);
  vertex(a/2, -b/2, c/2);
  vertex(-a/2, -b/2, c/2);
  endShape(CLOSE);
  
  fill(0,0,255);
  beginShape(QUAD);
  vertex(-a/2, -b/2, c/2);
  vertex(a/2, -b/2, c/2);
  vertex(a/2, b/2, c/2);
  vertex(-a/2, b/2, c/2);
  endShape(CLOSE);
  
  fill(0,0,128);
  beginShape(QUAD);
  vertex(-a/2, -b/2, -c/2);
  vertex(a/2, -b/2, -c/2);
  vertex(a/2, b/2, -c/2);
  vertex(-a/2, b/2, -c/2);
  endShape(CLOSE);
  
  stroke(255);
  noFill();
  for(int i = 0; i < engines.length; i++) {
    pushMatrix();
    translate((float)engines[i].position.x, (float)engines[i].position.y, (float)engines[i].position.z);
    box(0.05);
    popMatrix();
  }
  
  strokeWeight(4);
  stroke(255, 128, 0);
  for(int i = 0; i < body.forces.length; i++) {
    if(!body.forceActive[i]) continue;
    Vector3 start = body.forcePositions[i];
    Vector3 end = start.sub(body.forces[i]);
    beginShape(LINES);
    vertex((float)start.x, (float)start.y, (float)start.z);
    vertex((float)end.x, (float)end.y, (float)end.z);
    endShape();
  }
  
  fill(255);
  hint(DISABLE_DEPTH_MASK);
  perspective();
  camera();
  textFont(font);
  String[] viewModes = {"Orbit (body)", "Orbit (origin)", "FPV"};
  text(
    "CONTROLS:\n" + 
    "Mouse dragging: change view direction (orbit view)\n" +
    "Mousewheel: view closer/further (orbit view)\n" + 
    "V: cycle view - " + viewModes[viewMode] + "\n" + 
    "W,A,S,D: rotational RCS\n" +
    "H,N: forward/backward RCS\n" + 
    "I,J,K,L: up, down, left, right RCS\n" + 
    "T: Stabilization - " + (stabilizationEnabled ? "Enabled" : "Disabled") + "\n" +
    ",. TimeWarp x" + iterationsPerGUI, 20, 30);
  hint(ENABLE_DEPTH_MASK);
  
}

Vector3 getOrbitalVelocity(Vector3 position) {
  Vector3 direction = cross(new Vector3(0,0,1), position).unit();
  double speed = Math.sqrt(stdGravParam/position.mag());
  return direction.mult(speed);
}

void physicsUpdate() {
  stabilize(dt);
  body.forces = new Vector3[engines.length + 1];
  body.forceActive = new boolean[engines.length + 1];
  body.forcePositions = new Vector3[engines.length + 1];
  for(int i = 0; i < engines.length; i++) {
    engines[i].update(dt);
    body.forces[i] = engines[i].thrustVector;
    body.forceActive[i] = engines[i].active;
    body.forcePositions[i] = engines[i].position;
  }
  
  // Gravity
  Vector3 r = body.position;
  Vector3 gravity = r.unit().mult(-stdGravParam * body.mass / r.sqr());
  int index = body.forces.length - 1;
  body.forces[index] = body.basis.Transpose().mult(gravity);
  body.forceActive[index] = true;
  body.forcePositions[index] = new Vector3();
  
  body.update(dt);
}

void mouseDragged() {
  heading += -(mouseX - pmouseX) * 0.01;
  pitch += -(mouseY - pmouseY) * 0.01;
  pitch = constrain(pitch, -PI/2, PI/2);
}

void mouseWheel(MouseEvent event) {
  distance *= exp(event.getCount() * 0.1);
}

void manageInput() {
  // Translational
  engines[0].active = keys['w'];
  engines[1].active  = keys['a'];
  engines[2].active  = keys['s'];
  engines[3].active  = keys['d'];
  engines[4].active  = keys['n'];
  engines[5].active  = keys['s'];
  engines[6].active  = keys['d'];
  engines[7].active  = keys['w'];
  engines[8].active  = keys['a'];
  engines[9].active  = keys['h'];
  engines[10].active  = keys['q'];
  engines[11].active  = keys['e'];
  engines[12].active  = keys['l'];
  engines[13].active  = keys['e'];
  engines[14].active  = keys['q'];
  engines[15].active  = keys['j'];
  engines[16].active  = keys['i'];
  engines[17].active  = keys['k'];
}

boolean[] keys = new boolean[1024];
void keyPressed() {
  if(key < keys.length) keys[key] = true;
  if(key == 'v') viewMode = (viewMode + 1) % 3;
  if(key == 't') stabilizationEnabled = !stabilizationEnabled;
  if(key == '.') iterationsPerGUI = constrain(iterationsPerGUI + 1, 1, 10);
  if(key == ',') iterationsPerGUI = constrain(iterationsPerGUI - 1, 1, 10);
}

void keyReleased() {
  keys[key] = false;
}
