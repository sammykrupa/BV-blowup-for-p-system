%This is the MATLAB code to find a large T_N configuration in the set of
%2x2 matrices
%It is associated with the paper 'Finite time BV blow-up for Liu-admissible
%solutions to $p$-system via computer-assisted proof' by Sam G. Krupa


%Step 1: Setup the numerical solver to find a large T_N, and run the solver

exitflag=0;


%keep looping (with different random initial conditions for the solver) until the solver algorithm finds a viable solution
while exitflag~=1
    


%the epsilon determines the strictness in strict inequalities
epsilon_fixed=1;



%pick N\geq4 for the T_N configuration we are looking for
N=5;

%the kappa_i are random
%kappa=randi([2,10],N,1);

%to find the T_N we use in this paper, use
kappa=[19;37;67;64;73];

%define the optimization variables

deriv=optimvar('deriv',N,1,1,'UpperBound',-epsilon_fixed);

n=optimvar('n',1,2,N);
a=optimvar('a',2,1,N);

%define the optimization problem

prob=optimproblem;


%define the C_i rank-one matrices

C=optimexpr(2,2,N);
for i=1:N
C(:,:,i)=a(:,:,i)*n(:,:,i);
end



%we want the C_i to sum to zero (so sum along the third component of the C
%array)
prob.Constraints.C_sum_to_zero=sum(C(:,:,:),3)==zeros(2,2);



%we want the X(:,:,i) to be in T_N configuration
X=optimexpr(2,2,N);

for i=1:N
    X(:,:,i)=sum(C(:,:,1:(i-1)),3)+kappa(i)*C(:,:,i);
end



%we want function -p to be strictly decreasing and strictly convex
%(note the upper bound on the optimization variable deriv which is defined
%above)
for i=1:N
for j=[1:(i-1) (i+1):N]
constraint = strcat('FluxConstrConvexityCondition','_',num2str(i),'_',num2str(j));
prob.Constraints.(constraint)=X(2,2,j)-X(2,2,i)-deriv(i,1,1)*(X(1,1,j)-X(1,1,i))>=epsilon_fixed;
end
end

%we want the following symmetry-type condition on the T_N
for i=1:N
constraint = strcat('SymmetryCondition1','_',num2str(i));
prob.Constraints.(constraint)=X(2,1,i)+X(1,2,i)==0;
end


%initialize the solver with random initial points

%data.a = randi([-1000,1000],2,1,N);
%data.n=randi([-1000,1000],1,2,N);
%data.deriv=randi([-1000,-2],N,1,1);

%to find the T_N we use in this paper, use
data.n(:,:,1)=[518 729];
data.n(:,:,2)=[714 -734];
data.n(:,:,3)=[727 189];
data.n(:,:,4)=[181 836];
data.n(:,:,5)=[-889 410];

data.a(:,:,1)=[345;-915];
data.a(:,:,2)=[-401;-25];
data.a(:,:,3)=[-312;-849];
data.a(:,:,4)=[-328;-755];
data.a(:,:,5)=[-715;-247];

data.deriv=[-503;-635;-3;-46;-488];



%solve!
disp('Start the solver')
[sol,fval,exitflag,output] = solve(prob,data,'Options',optimoptions(prob,'MaxFunctionEvaluations',200000,'MaxIterations',200,'ConstraintTolerance',.00000000000000001,'Display','iter'));
end %suitable solution found, end the while loop



%save the original solution (before we start modifying it). Just in case
%we want to examine it later.
original_sol=sol;

%book keeping
sol.kappa=kappa;

%convert all numeric values (from the numerical optimization) to symbolic expressions
%from now on, we only work symbolically (not numerically)
sol.deriv=sym(sol.deriv,'f');
sol.a=sym(sol.a,'f');
sol.n=sym(sol.n,'f');
sol.P(:,:)=sym(zeros(2,2),'f'); %just for book-keeping, make a note that P=0;
sol.kappa=sym(kappa,'f');

%Step 2: calculate the X_i and check the algebraic constraints.

%now, calculate the X_i (making up the T_N) and calculate the rank-one matrices C_i. 
X=sym(zeros(2,2,N));
C=sym(zeros(2,2,N));
for i=1:N
C(:,:,i)=sol.a(:,:,i)*sol.n(:,:,i);
end
for i=1:N
X(:,:,i)=sum(C(:,:,1:(i-1)),3)+sol.kappa(i)*C(:,:,i);
end

%FORCE the symmetry-type condition to be EXACTLY true
for i=1:N
X(2,1,i)=-X(1,2,i);
end

%check the algebraic conditions on the T_N to make sure they are satisfied


%we want function -p to be strictly decreasing and strictly convex
%we will approximate symbolic values to 32 significant digits using variable-precision arithmetic (arbitrary-precision floating-point numbers)
%see https://www.mathworks.com/help/symbolic/vpa.html
disp('Check the FluxConstrConvexityCondition -- these must all be POSITIVE');

for i=1:N
for j=[1:(i-1) (i+1):N]
vpa(X(2,2,j)-X(2,2,i)-sol.deriv(i,1,1)*(X(1,1,j)-X(1,1,i)))
end
end

disp('Check the values of the deriv variable -- these must all be NEGATIVE');
vpa(sol.deriv)

%Now, we check that the T_N contains no rank-one connections
disp('Check that the T_N contains no rank-one connections...these rank values values should all be 2')
for i=1:N
for j=[1:(i-1) (i+1):N]
%calculate rank
rank(X(:,:,i)-X(:,:,j))
end
end




%Step 3: check the various orderings of the T_N

%for each ordering, we compute the paramtrization (P,C_i,kappa_i) (following
%Székelyhidi, László, Jr.: Rank-one convex hulls in R^{2×2}. Calc.Var.PDE 22(3),253–281(2005))
%see also Förster, Clemens and Székelyhidi, László, Jr.: T5-configurations and non-rigid sets of matrices. Calc. Var. Partial Differential Equations 57(1), Art. 19 (2018)
%We use the orderings which we determined by trial and error for the particular large T_5 used in this paper:
SavedOrderings=[1:N;4 3 1 5 2;4 2 3 1 5]; %one ordering for each row
%initialize the symbolic matrices we will use to store our
%parameterizations
sol.P=sym(zeros(2,2,3));
sol.C=sym(zeros(2,2,N,3));
sol.kappa=sym(zeros(N,1,3));

for q=1:3 %loop through the three orderings

%symbolically find the scalar mu and the vector lambda 

syms mu;

    for i=1:N
    for j=1:N
        temp=X(:,:,SavedOrderings(q,i))-X(:,:,SavedOrderings(q,j));
        if i>j
            %compute determinant of temp matrix (multiplied by mu)
            Amu(i,j)=mu*(det(temp));
        else
            %compute determinant of temp matrix
            Amu(i,j)=det(temp);
        end
    end
    end

%solve for mu symbolically using root function 
%see https://www.mathworks.com/help/symbolic/sym.root.html

zerosVal=root(det(Amu),mu);
muValue=zerosVal(4); %for the given orderings in SavedOrdering, the fourth zero always give an acceptable value of mu (>1)
%substitute this value of mu into the matrix Amu, then calculate a basis
%for the null space of Amu
basis=null(subs(Amu,mu,muValue)); %for the given orderings in SavedOrdering, this basis is one vector, each component of which is positive. Thus, we can take basis=lambda;
lambda=basis; %the basis is in fact the lambda vector from the algebraic criterion for a T_N
 

%Step 4: For each ordering, compute the parameterization of the T_N. Then, check the
%linear independence of the rank-one directions given by the C_i.

%now, from the algebraic criterion, we can compute the parameterization
%(P,C_i,kappa_i) for this ordering of the T_N

 C=sym(zeros(2,2,N));
 kappa=sym(zeros(N,1));
 tempP=sym(zeros(2,2,N)); %for each ordering of the T_N, this will store the P_i associated with the vertices of the inner `rank-one' loop


    %compute the P_i
    for i=1:N
            %compute the scalar pre-factor
            denominator=sum(lambda(i:N,1))+muValue*sum(lambda(1:(i-1),1));
            %compute the matrix term
            matrixSum=sym(zeros(2,2));
            for k=1:N
                if k<i
                    matrixSum=matrixSum+muValue*lambda(k,1)*X(:,:,SavedOrderings(q,k));
                else
                    matrixSum=matrixSum+lambda(k,1)*X(:,:,SavedOrderings(q,k));
                end
            end
   tempP(:,:,i)=(1/denominator)*matrixSum;

         if i==1
           sol.P(:,:,q)=tempP(:,:,i); %save the initial P values for the first, second and third orderings
         end
         
         if i>1
             tempC=tempP(:,:,i)-tempP(:,:,i-1); %this gives the C_{i-1}
             tempkappaC=X(:,:,SavedOrderings(q,i-1))-tempP(:,:,i-1); %this gives kappa_i C_{i-1}
             kappa(i-1,1)=tempkappaC(1,1)/tempC(1,1); %assuming tempC(1,1) is not zero
             C(:,:,i-1)=tempC; %save the C_i
         end

    end
    %now, collect the info from the last of the P_i (i=N) after the for
    %loop above has finished
             tempC=tempP(:,:,1)-tempP(:,:,N); %this gives the C_N
             tempkappaC=X(:,:,SavedOrderings(q,N))-tempP(:,:,N); %this gives kappa_N C_N
             kappa(N,1)=tempkappaC(1,1)/tempC(1,1);
             C(:,:,N)=tempC; %save the C_N

              %save all of our results from determining the
              %parameterization corresponding to this particular ordering
              %of the T_N (ordering #q)
             sol.C(:,:,:,q)=C;
             sol.kappa(:,:,q)=kappa;
end



%check the linear independence of the rank-one directions as required by
%the definition of `large' T_5

   disp('Check the linear independence of rank-one directions -- these determinants must all be NON-ZERO');

for i=1:N

   calculateRank=sym(zeros(3,4));
  
   calculateRank=[reshape(sol.C(:,:,i,1),1,4);reshape(sol.C(:,:,i,2),1,4);reshape(sol.C(:,:,i,3),1,4)];
   %determine the exact rank of the symbolic matrix calculateRank by
   %showing a non-zero 3x3 minor (using symbolic computation of
   %determinant), which shows rank = 3
   vpa(det(calculateRank(1:3,[1 3 4])))
end
