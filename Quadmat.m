% Quadcopter model and simulation:
%__________________________________________________________________________


%__________________________________________________________________________
% Simulation time: 
    dt=0.01;
    ts=0:dt:10; % time sequence
    N=numel(ts);  % Nr of elemements in the sequence
%__________________________________________________________________________


%__________________________________________________________________________
% Physical properties of the Quadcopter and its enviroment:
    
    g=9.82;                                         % gravity
    m=0.2;                                          % mass

    Jxx=1;                                          % Inertia around x axis
    Jyy=1;                                          % Inertia around y axis
    Jzz=1;                                          % Inertia around z axis
    J=diag([Jxx,Jyy,Jzz]);                          % inertia matrix
    
    b=1e-7;                                         % Lift coefficient of engine
    k=3e-6;                                          % Air resistance coefficient
    d=0.25;                                         % distance to engines from C.G
    
    Kd=0.25;                                        % drag force coefficient
%__________________________________________________________________________

%__________________________________________________________________________
% Init empty arrays to save state information
    pos= zeros(3,N);        % position
    vel=zeros(3,N);         % velocity
    angls=zeros(3,N);       % angles
    angvel=zeros(3,N);      % angle velocities
    in=zeros(4,N);                % input
%__________________________________________________________________________


%__________________________________________________________________________
% Initial system state

    xyz=zeros(3,1);                 % starting postition vector
    xyzdot=zeros(3,1);              % Initial velocity vector
    
    Ax=0;                           % Initial angle around x axis
    Ay=0;                           % Initial angle around y axis
    Az=0;                           % Initial angle around z axis
    
    As=[Ax;Ay;Az];
    
    Asdot=zeros(3,1);               % Initial angular velocities
%__________________________________________________________________________

%__________________________________________________________________________
%INPUTS:
    % Test inputs:
        %  Angular Velocity of engines:
            w1=0;
            w2=0;
            w3=0;
            w4=0;
            Ave=[w1;w2;w3;w4];
            Inputs = Ave.^2;
%__________________________________________________________________________

index=0;

dc=zeros(3,1); %derivative part of controller
ic=zeros(3,1); % integral part of controller

kd=0.5;
kp=0.01;
ki=0.1;
for t=ts
  index=index+1;
  
  % Controller:.....
  % In : Asdot, + parameters
  % Out: new controll inputs + parameters
  %___________________________________________________
  
  % PID
  
  % To prevent wind up
  if (max(abs(ic)))> 0.01
      ic(:)=0;
  end
  
  % Total thrust PID
   TotalThrust= m*g/(b*(cos(dc(1))*cos(dc(2))));
   
   error= kd*Asdot + kp * dc - ki *ic;
   Inputs(1)=TotalThrust/4 -(2*k *error(1)*Jxx+ error(3) * Jzz * b *d ) / (4*b*k*d);
   Inputs(2)=TotalThrust/4 + error(3) * Jzz/(4*k) - (error(2) *Jyy) / (2*b*d);
   Inputs(3)=TotalThrust/4 -(-2 *b *error(1) *Jxx + error(3) *Jzz *b *d )/(4* k *b *d);
   Inputs(4)=TotalThrust/4 + error(3) * Jzz /(4*k) + (error(2)*Jyy)/(2 *b *d);
   
   dc= dc + dt.* Asdot;
   ic= ic+ dt.*dc;
   %__________________________________________________
   
  % W -omega- angular velocity wrt the individual axes
  % WR- Transformation matrix: 
  WR= [1        0                -sin(Ay)         ;
       0      cos(Ax)          cos(Ay)*sin(Ax)    ;
       0     -sin(Ax)          cos(Ay)*cos(Ax)    ];
            
  W = WR * Asdot;
  %___________________________________________________
  
  % engine constants relating lift and torque to angular velocity
            C= [-b  -b  -b  -b;
                 0 -d*b  0  d*b;
                d*b  0 -d*b  0;
                 k  -k   k  -k];
  
      Tt=C*Inputs; % Total thrust and torques [T; tx; ty; tz]
  
  % Acceleration of inertial ref frame
      Fd=Kd * xyzdot; % Drag force

      % Rotation matrices:
                Rx=[1   0       0     ;
                    0 cos(Ax) -sin(Ax);
                    0 sin(Ax)  cos(Ax)];

                Ry=[cos(Ay)  0 sin(Ay);
                    0        1       0;
                    -sin(Ay) 0 cos(Ay)];

                Rz=[cos(Az) -sin(Az) 0;
                    sin(Az)  cos(Az) 0;
                    0           0    1];

                % Total Rotation matrix R (XYZ convention)
                R=Rz*Ry*Rx;
      
      
      acceleration= [0;0;g]- R*[0;0;Tt(1)]/m -Fd/m; 
  %___________________________________________________
  % Angular acceleration in body frame 
         dW=J^(-1)*(-cross(W, J*W))+Tt(2:4);
  
  % Update states
   W = W+ dt*dW;
   Asdot=WR^(-1)*W;
   As = As + dt* Asdot;
   xyz= xyz + dt*xyzdot;
   
   % Save data
    pos(:,index)= xyz;              % position
    vel(:,index)=xyzdot;            % velocity
    angls(:,index)=As;              % angles
    angvel(:,index)=Asdot;          % angle velocities
    in(:,index)=Inputs;             % input  
    
end    