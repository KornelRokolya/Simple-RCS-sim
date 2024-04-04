boolean stabilizationEnabled = false;
Vector3 targetOmega = new Vector3();
Matrix3x3 targetOrientation = Identity3x3();
PID xPID = new PID(20, 0, 1, 0.5);
PID yPID = new PID(20, 0, 1, 0.5);
PID zPID = new PID(20, 0, 1, 0.5);

void stabilize(double dt) {
  for(int i = 0; i < engines.length; i++) engines[i].throttle = 1;
  for(int i = 0; i < keys.length; i++) {
    if(keys[i]) {
      if(keys['w'] || keys['a'] || keys['s'] || keys['d'] || keys['q'] || keys['e']) {
        targetOrientation = body.basis.Copy();
      }
      return;
    }
  }
  if(!stabilizationEnabled) return;
  
  Matrix3x3 rotation = targetOrientation.mult(body.basis.Transpose());
  double angle = Math.acos(constrain((rotation.Trace() - 1) / 2, -1, 1));
  targetOmega = body.basis.Transpose().mult(
      new Vector3(
          rotation.a21 - rotation.a12,
          rotation.a02 - rotation.a20,
          rotation.a10 - rotation.a01));
  if(angle != 0) targetOmega = targetOmega.mult(angle / Math.sin(angle));
  if(Math.abs(angle) < 0.05) targetOmega = targetOmega.mult(angle / Math.sin(angle) * 0.1);
  //println(angle);
  
  Vector3 omegaLocal = body.basis.Transpose().mult(body.angularVelocity);
  xPID.update(omegaLocal.x, targetOmega.x, dt);
  yPID.update(omegaLocal.y, targetOmega.y, dt);
  zPID.update(omegaLocal.z, targetOmega.z, dt);
  
  // Set all engines to inactive
  for(int i = 0; i < engines.length; i++) {
    engines[i].active = false;
  }
  
  // X output (roll)
  double X = xPID.output;
  if(X >= 0) {
    engines[10].active = true;
    engines[10].throttle = X;
    engines[14].active = true;
    engines[14].throttle = X;
  }
  else {
    engines[11].active = true;
    engines[11].throttle = -X;
    engines[13].active = true;
    engines[13].throttle = -X;
  }
  
  // Y output (pitch)
  double Y = yPID.output;
  if(Y >= 0) {
    engines[2].active = true;
    engines[2].throttle = Y;
    engines[5].active = true;
    engines[5].throttle = Y;
  }
  else {
    engines[0].active = true;
    engines[0].throttle = -Y;
    engines[7].active = true;
    engines[7].throttle = -Y;
  }
  
  // Z output (heading)
  double Z = zPID.output;
  if(Z >= 0) {
    engines[1].active = true;
    engines[1].throttle = Z;
    engines[8].active = true;
    engines[8].throttle = Z;
  }
  else {
    engines[3].active = true;
    engines[3].throttle = -Z;
    engines[6].active = true;
    engines[6].throttle = -Z;
  }
}

class PID {
    public double value, valuePrev;
    public double reference;
    public double error, errorPrev;
    public double errorI;
    public double errorD;
    public double P, I, D, TD2;
    public double Pmultiplier = 1;
    public double Imultiplier = 1;
    public double Dmultiplier = 1;
    public double output;

    public PID(double P, double I, double D, double TD2) {
        this.P = P;
        this.I = I;
        this.D = D;
        this.TD2 = TD2;
    }

    public void Reset() {
      errorI = 0;
    }

    public void update(double value, double ref, double dt) {
        valuePrev = value;
        this.value = value;
        this.reference = ref;
        errorPrev = error;
        error = reference - value;
        errorI += error * dt;
        if(dt != 0) errorD += ((error - errorPrev) / dt - errorD) * dt / (TD2 + dt);
        output = P * Pmultiplier * error + I * Imultiplier * errorI + D * Dmultiplier * errorD;
    }
}
