model quad
// Quadcopter model and simulation:
//__________________________________________________________________________
  
 constant Integer N=1001;
 constant Real dt=0.01;
          Real ts=0;
//__________________________________________________________________________
// Physical properties of the Quadcopter and its enviroment:
    
   parameter Real g=9.82;                                         // gravity
   parameter Real m=0.2;                                          // mass

   parameter Real Jxx=1;                                          // Inertia around x axis
   parameter Real Jyy=1;                                          // Inertia around y axis
   parameter Real Jzz=1;                                          // Inertia around z axis
   parameter Real J [3,3]=[{Jxx,0,0}, {0,Jyy,0},{0,0,Jzz}];                          // inertia matrix
   
   parameter Real b=1e-7;                                         // Lift coefficient of engine
   parameter Real k=3e-6;                                          // Air resistance coefficient
   parameter Real d=0.25;                                         // distance to engines from C.G
    
   parameter Real Kd=0.25;                                        // drag force coefficient
//__________________________________________________________________________


//__________________________________________________________________________
// Init empty arrays to save state information
   parameter Real pos[3,N]= zeros(3,N);        // position
   parameter Real vel[3,N]=zeros(3,N);         // velocity
   parameter Real angls[3,N]=zeros(3,N);       // angles
   parameter Real angvel[3,N]=zeros(3,N);      // angle velocities
   parameter Real inp [4,N]=zeros(4,N);                // input
//________________________________________________________________  __________

//__________________________________________________________________________
// Initial system state

   parameter Real xyz[3]=zeros(3);                 // starting postition vector
   parameter Real xyzdot[3]=zeros(3);              // Initial velocity vector
    
   parameter Real Ax=0;                           // Initial angle around x axis
   parameter Real Ay=0;                           // Initial angle around y axis
   parameter Real Az=0;                           // Initial angle around z axis
    
   parameter Real As [3]={Ax,Ay,Az};
    
   parameter Real Asdot[3]=zeros(3);               // Initial angular velocities
//__________________________________________________________________________

//__________________________________________________________________________
//INPUTS:
    // Test inputs:
        //  Angular Velocity of engines:
           parameter Real w1=0;
           parameter Real w2=0;
           parameter Real w3=0;
           parameter Real w4=0;
           parameter Real Ave [4]={w1,w2,w3,w4};
           parameter Real Inputs [4]=Ave;
           //Real Inputs = Ave.^2;
//__________________________________________________________________________

// Init contol parameters 
  parameter Real index=0;
  
  parameter Real dc[3]=zeros(3); //derivative part of controller
  parameter Real ic[3]=zeros(3); // integral part of controller
  
  parameter Real kd=0.5;
  parameter Real kp=0.01;
  parameter Real ki=0.1;
  parameter Real PidThrust;
  parameter Real error[3];
  parameter Real WR[3,3];
  parameter Real W[3];
  parameter Real C[4,4];
  parameter Real Tt[4];
  parameter Real Fd[3];
  parameter Real Rx[3,3];
  parameter Real Ry[3,3];
  parameter Real Rz[3,3];
  parameter Real R [3,3];
  parameter Real dW[3];
  
  equation
  for index in 1:N loop
     ts= dt*index;
    
      // Controller:.....
      // In : Asdot, + parameters
      // Out: new controll inputs + parameters
      //___________________________________________________
      
      // PID
      // To prevent wind up
      if max(abs(ic))>0.01 then
        ic=zeros(3);
      end if;

      // Total thrust PID
      PidThrust= m*g/(b*(cos(dc[1])*cos(dc[2])));
      error= kd*Asdot + kp * dc - ki *ic;
      Inputs[1]=PidThrust/4 -(2*k *error[1]*Jxx+ error[3] * Jzz * b *d ) / (4*b*k*d);
      Inputs[2]=PidThrust/4 + error[3] * Jzz/(4*k) - (error[2] *Jyy) / (2*b*d);
      Inputs[3]=PidThrust/4 -(-2 *b *error[1] *Jxx + error[3] *Jzz *b *d )/(4* k *b *d);
      Inputs[4]=PidThrust/4 + error[3] * Jzz /(4*k) + (error[2]*Jyy)/(2 *b *d);
      
      dc= dc + dt.* Asdot;
      ic= ic+ dt.*dc;
      
      //__________________________________________________
   
      // W -omega- angular velocity wrt the individual axes
      // WR- Transformation matrix: 
             WR= [{1    ,    0          ,      -sin(Ay)},         
                  {0    ,  cos(Ax)      ,    cos(Ay)*sin(Ax)},    
                  {0    , -sin(Ax)      ,    cos(Ay)*cos(Ax)}];
            
             W = WR * Asdot;
             
       //___________________________________________________
  
       // engine constants relating lift and torque to angular velocity
             C= [{-b  , -b  , -b  , -b },
                 { 0  ,-d*b ,  0  , d*b},
                 {d*b ,  0  ,-d*b ,  0 },
                 { k  , -k  ,  k  , -k }];
                 
             Tt=C*Inputs; // Total thrust and torques [T; tx; ty; tz]
             
       // Acceleration of inertial ref frame
             Fd=Kd * xyzdot; // Drag force  
             
             // Rotation matrices:
                Rx=[{1 ,  0    ,    0    },
                    {0 ,cos(Ax), -sin(Ax)},
                    {0 ,sin(Ax),  cos(Ax)}];

                Ry=[{cos(Ay) ,  0, sin(Ay)},
                    {0       ,  1,    0   },
                    {-sin(Ay),  0, cos(Ay)}];

                Rz=[{cos(Az), -sin(Az), 0},
                    {sin(Az),  cos(Az), 0},
                    {0      ,     0   , 1}]; 
              // Total Rotation matrix R (ZYX convention)
                R=Rz*Ry*Rx;   
                
              //___________________________________________________
              // Angular acceleration in body frame
              
              //W3=cross(W3, W3);
                dW= Modelica.Math.Matrices.inv(J)*(-cross(W, J*W))+Tt[2:4];  
                
              // Update states
                W = W+ dt*dW;
                Asdot=Modelica.Math.Matrices.inv(WR)*W;
                As = As + dt* Asdot;
                xyz= xyz + dt*xyzdot;  
                
                // Save data
                pos[:,index]= xyz;              // position
                vel[:,index]=xyzdot;            // velocity
                angls[:,index]=As;              // angles
                angvel[:,index]=Asdot;          // angle velocities
                inp[:,index]=Inputs;            // input  
  end for;
end quad;