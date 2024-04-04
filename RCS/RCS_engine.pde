RCS_engine[] engines;

public class RCS_engine {
  public Vector3 position;
  public Vector3 direction;
  public Vector3 thrustVector;
  public double maxThrust = 1;
  public double throttle = 1;
  public boolean active = false;
  public boolean pulseModeEnabled = true;
  public double pulseFrequency = 10;  
  double pulsePhase = 0;
  
  public RCS_engine(Vector3 pos, Vector3 dir) {
    this.position = pos;
    this.direction = dir;
    thrustVector = new Vector3();
  }
  
  public RCS_engine(double dx, double dy, double dz, double x, double y, double z) {
    this.position = new Vector3(x,y,z);
    this.direction = new Vector3(dx,dy,dz).unit();
    thrustVector = new Vector3();
  }
  
  public void update(double dt) {
    double thr = maxThrust;
    if(!active) thr = 0;
    throttle = constrain((float)throttle, 0, 1);
    
    if(pulseModeEnabled) {
      pulsePhase = (pulsePhase + pulseFrequency * dt) % 1;
      if(pulsePhase > throttle) thr = 0;
    }
    else {
      pulsePhase = 0;
      thr *= throttle;
    }
    
    thrustVector = direction.mult(thr);
  }
}
