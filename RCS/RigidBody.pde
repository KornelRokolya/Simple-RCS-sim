class RigidBody {
  public double mass;
  public Matrix3x3 Theta0;
  public Matrix3x3 basis;
  public Vector3 position;
  public Vector3 velocity;
  public Vector3 angularMoment;
  public Vector3 angularVelocity;
  public Vector3 acceleration;
  
  public Vector3[] forces;
  public Vector3[] forcePositions;
  public boolean[] forceActive;
  public Vector3[] torque;
  public boolean[] torqueActive;
  
  public void update(double dt) {
    Vector3 netForce = new Vector3();
    Vector3 netTorque = new Vector3();
    
    for(int i = 0; i < forces.length; i++) {
      if(!forceActive[i]) continue;
      netForce.addTo(forces[i]);
      netTorque = netTorque.add(cross(forcePositions[i], forces[i]));
    }
    
    for(int i = 0; i < torque.length; i++) {
      if(!torqueActive[i]) continue;
      netTorque.addTo(torque[i]);
    }
    
    netForce = basis.mult(netForce);
    netTorque = basis.mult(netTorque);
    
    
    // Positional update
    Vector3 acc = netForce.div(mass);
    position = lin(position, 1, velocity, dt, acc, dt*dt/2);
    velocity = lin(velocity, 1, acc, dt);
    acceleration = acc;
    
    // Rotational update
    Matrix3x3 Theta = basis.mult(Theta0.mult(basis.Transpose()));
    Vector3 omega = Theta.Inverse().mult(angularMoment);
    Vector3 n = omega.unit();
    double angle = omega.mag() * dt;
    Matrix3x3 dR = GetRotationMatrix_unitAxisAngle(n, angle);
    basis = dR.mult(basis);
    
    angularMoment = lin(angularMoment, 1, netTorque, dt);
    
    if(position.sqr() <= R*R) {
      position = position.unit().mult(R);
      velocity = new Vector3();
      angularMoment = new Vector3();
    }
    
    // Calculate omega
    Theta = basis.mult(Theta0.mult(basis.Transpose()));
    angularVelocity = Theta.Inverse().mult(angularMoment);
  }
}
