%This illustrates the overlapping grid method for solving diffusion
%equations

clc
clear all 
%Assuming number of grid points in Fine grids are J and Number of grid
%points in coarse grids is K 

fine_mesh = 0.01 ; %meshsize for the finer grid 
coarse_mesh = 0.1 ; % meshsize for the coarse grid 

%Starting of the overlap 
overlap_cell =  50; %start of the coarse grid for overlapping 

J = 60 ; %Number of Control Volumnes in Fine grid  
K = (1 - (overlap_cell * fine_mesh))/coarse_mesh  %Number of Control Volumnes in Coarse grid 

%Domain Description 
x_s = 0 ; 
x_f = 1 ; 

%Indentify the location of the interpolation point in both coarse and fine
%grid 
if (fine_mesh*J >(x_f- coarse_mesh*K) && (fine_mesh*J < x_f && coarse_mesh*K < x_f) )
    
    %finding neighbouring points in coarse grid 
    fine_grid_ghost_point = fine_mesh * J + fine_mesh*0.5 
    coarse_grid_cell = floor(((1- fine_grid_ghost_point)/coarse_mesh))
    left_over_dist = (1-fine_grid_ghost_point) - coarse_grid_cell * coarse_mesh
    
    %For interpolant ahead of half point 
    if(left_over_dist > (0.5*coarse_mesh))
        coarse_grid_point_l = coarse_grid_cell+2    
        coarse_grid_point_r = coarse_grid_cell+1  
        
    elseif(left_over_dist ==0)
        coarse_grid_cell = floor(((1- fine_grid_ghost_point)/coarse_mesh))-1;
        coarse_grid_point_l = coarse_grid_cell + 2 
        coarse_grid_point_r = coarse_grid_cell  + 1 
        
    %for interpolant point lagging cell center   
    elseif (left_over_dist<(0.5*coarse_mesh)) 
        coarse_grid_point_l = coarse_grid_cell+1  
        coarse_grid_point_r = coarse_grid_cell     
        
    %for interpolant ahead at cell center
    elseif (left_over_dist==(0.5*coarse_mesh))
        coarse_grid_point_l = coarse_grid_cell 
        coarse_grid_point_r = coarse_grid_cell+1   

    %if interpolant lies at the wall of cell  

    end

    
    %Finding Neighboring points in Fine grids
    coarse_grid_ghost_point = overlap_cell * fine_mesh - (0.5*coarse_mesh)
    fine_grid_cell = floor(coarse_grid_ghost_point/fine_mesh) 
    left_over_dist =  coarse_grid_ghost_point - fine_grid_cell*fine_mesh  
    
    if(left_over_dist > 0.5*fine_mesh)
        fine_grid_point_l = fine_grid_cell+1 
        fine_grid_point_r = fine_grid_cell+2  
        
    elseif (left_over_dist== (0.5*fine_mesh))
        fine_grid_point_l = fine_grid_cell+1
        fine_grid_point_r = fine_grid_cell+2
        
    elseif(left_over_dist==0)
        fine_grid_cell = floor(coarse_grid_ghost_point/fine_mesh) -1;
        fine_grid_point_l = fine_grid_cell+1 ;
        fine_grid_point_r = fine_grid_cell+2 ;
        
    elseif(left_over_dist< 0.5 * fine_mesh)
        fine_grid_point_l = fine_grid_cell
        fine_grid_point_r = fine_grid_cell+1 
    end
    
else 
    display('Something Went wrong,there is no overlapping ')
    quit cancel
end

% If we have reached this point , that means there is overlapping 
% We have to find the distances of the interpolant with neighbouring points 
% For fine grid ghost point - interpolant for the fine grids 

x2 = (1-fine_grid_ghost_point)-(coarse_grid_cell+0.5)*coarse_mesh 
x1 = coarse_mesh-x2
beta_fine_l = x2/coarse_mesh 
beta_fine_r = x1/coarse_mesh  

%for the coarse grid points - interpolant for the coarse grids 
x1 = coarse_grid_ghost_point - (fine_grid_cell+0.5)*fine_mesh 
x2 = fine_mesh-x1  
beta_coarse_l = x2/fine_mesh 
beta_coarse_r = x1/fine_mesh

N = J + K ;
A = zeros(N,N) ; 

%Matrix has to be assembled 

%Matrix Coefficient 
alpha_l = -1 / fine_mesh ; 
alpha_c = 2 / fine_mesh ; 
alpha_r = -1 / fine_mesh ; 

beta_l = -1 / coarse_mesh ; 
beta_c = 2 / coarse_mesh ; 
beta_r = -1 / coarse_mesh ;  

A(1,1) = 3/fine_mesh ;
A(1,2) = -1/fine_mesh ; 
A(N,N) = 3/coarse_mesh;
A(N,N-1) = -1/coarse_mesh;
% 
for i = 2 : J-1
    
    %Assembling the matrix for the fine part except Ghost point and
    %Boundary Points 
    A(i,i-1) = alpha_l ;  
    A(i,i) = alpha_c ;
    A(i,i+1) = alpha_r ; 
end

for i = J+2 : N-1 

    %Assembling the matric for coarse part except the Ghost Point and
    %Boundary Poitns
    A(i,i-1) = beta_l ; 
    A(i,i)   = beta_c ; 
    A(i,i+1) = beta_r ; 
    
end

%Matrix entry for the ghost points on both grids 
A(J,J) = 2/fine_mesh ;
A(J,J-1)= -1/fine_mesh; 
A(J,K-coarse_grid_point_l+J+1) = -beta_fine_l / fine_mesh;
A(J,K-coarse_grid_point_r+J+1) = -beta_fine_r / fine_mesh; 

A(J+1,J+1) = 2/coarse_mesh ;
A(J+1,J+2) = -1/coarse_mesh;
A(J+1,fine_grid_point_l) = -beta_coarse_l/coarse_mesh ;
A(J+1,fine_grid_point_r) = -beta_coarse_r/coarse_mesh ; 

%Boundary Conditions 
a = 0.1 ; 
b = 1 ; 

bd = zeros(N,1) ;
gamma_0 = 2*a/fine_mesh ; 
gamma_1 = 2*b/coarse_mesh;
bd(1,1) = gamma_0 ; 
bd(N,1) = gamma_1 ; 

q = zeros(N,1) ; 

%Assembling source term for fine mesh 
for i = 1 : J 
    q(i) = fine_mesh* sin(pi * (i-0.5)*fine_mesh )  ; 
end

%Assembling source term for the coarse mesh 
for i = 1 : K 
    q(N-i+1) = coarse_mesh*  sin((x_f - (i-0.5)*coarse_mesh)*pi) ;  
end

r = bd + q ; 
u = A\r

fine_x = zeros(J,1) ; 
for i = 1 : J 
    fine_x(i) = (i-0.5) * fine_mesh ; 
end

coarse_x = zeros(K,1) ;
for i = 1 : K 
    coarse_x(i) = (x_f - (K-i+0.5)*coarse_mesh); 
end
fine_x 
coarse_x

sol_fine = zeros(J,1);
for i = 1 : J 
    sol_fine(i) = u(i) ;
end
sol_coarse = zeros(K,1) ;
for i = 1 : K 
    sol_coarse(i) = u(J+i) ; 
end
sol_fine 
sol_coarse

plot(fine_x,sol_fine,'color','r');
hold on ;
plot(coarse_x,sol_coarse,'color','b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TEST PROGRAM %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc
% clear all
% 
% N = 10 ; 
% h_e = 1/N ; 
% B = zeros(N,N) ; 
% 
% a = 0.1 ; 
% b = 1 ; 
% 
% for i = 2 : N-1 
%     B(i,i-1) = -1/h_e ; 
%     B(i,i) = 2/h_e ;
%     B(i,i+1) = -1/h_e ; 
% end
% 
% B(1,1) = 3/h_e ;
% B(1,2) = -1/h_e ; 
% B(N,N) = 3/h_e ; 
% B(N,N-1) = -1/h_e ; 
% 
% bd = zeros(N,1); 
% bd(1,1) = 2*a/h_e ; 
% bd(N,1) = 2*b/h_e ; 
% 
% q = zeros(N,1) ; 
% for i = 1 : N 
%     q(i) = h_e * sin((i-0.5)*h_e *pi) ; 
% end
% 
% r = bd + q ; 
% 
% u = B\r   
% i = 1:1:N; 
% plot(i,u)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TEST PROGRAM %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    


