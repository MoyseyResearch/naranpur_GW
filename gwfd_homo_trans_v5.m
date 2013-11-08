

%calculate parameters:
T = K*b;%H(:,ti);   %calculate the transmissivity
S = S*b; %Ss*b + Sy;%calculate storativity (currently fixed as parameter)

%% update times for error control

criterion = K*b/S*dt/dx^2;
disp(['Criterion with dt=',num2str(dt),' is ',num2str(criterion)])
if criterion > err_tol; 
    dt = err_tol*dx^2*S/K/b;
    disp(['Time step updated for improved performance, dt=',num2str(dt)])
end

t = tmin:dt:tmax;
Nt = length(t);


%% Actual code for simulation
p = S*dx^2/dt./T;        %evalute coef p for current H

%total number of cells
Ncells = Nx*Ny;

%high estimate of non-zero cells:
max_nonzeros = 5*Ncells;

%initialize matrices:
A = spalloc(Ncells,Ncells,max_nonzeros);  %initialize sparse matrix
rhs = zeros(Ncells,1);
H = zeros(Nx*Ny,Nt);    %initialize solution matrix
H(:,1) = Ho;

%generate the A matrix:
for i = 1:Ncells
    %get index for neighbor cells
    neighbors = [-Ny Ny -1 1] + i;
    %remove values outside of grid domain
    neighbors(neighbors<1) = [];
    neighbors(neighbors>Ncells) = [];
    
    if cells(i) == 1;  %this is an interior cell
        A(i,neighbors) = 1;
        A(i,i) = -(4+p);
    elseif cells(i) == 2;   %constant head cell
        A(i,i) = 1;
    elseif cells(i) == 3; %constant flux boundary
        for j=1:length(neighbors);
            if cells(neighbors(j))==1;
                A(i,neighbors(j)) = 1;
                A(i,i) = -1;
            end
        end
    end
end

%run through time steps
for ti = 1:Nt;
    rhs = -p*H(:,ti);       %set right hand side head from previous time step
    rhs = rhs-R*dx^2/T *(0.5*sin(ti/53*2*pi)+0.5); %add recharge rate (since used here for seasonality)
    rhs(pump_cell) = rhs(pump_cell)-Qwell(:)/T.*(well_elev<H(pump(:,1),ti))*(t(ti)>0.3*max(t)); %add wells  -- (t(ti)>0.3*max(t)) is just for testing to turn the wells on 1/3 of the way the into the simulation
    rhs(cells==2) = bH(cells==2);   %update head boundaries (overwrite rhs)
    rhs(cells==3) = -bq(cells==3)*dx; %update flux boundaries (overwrite rhs)
    H(:,ti+1) = bicg(A,rhs,1e-4,50);    %solve with conjugate gradients
    disp(['Done time: ',num2str(ti),' of ',num2str(Nt)])    %remove display line
end
