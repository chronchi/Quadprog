#Similar function to quadprog routine in MATLA. Here f(x) = (1/2)*x'*H*x + f'*x,
#where H is a matrix, f a vector and we have equality as well as
#boundary restrictions such as Ax <= b, (Aeq)*x = beq, lb <= x <= ub.

#This routine uses as basis de Ipopt Solver, speeding up the optimization.
#but one can use other solvers, just change the usesolver as input of the
#funtion

using JuMP, Ipopt

function quadprog(H,f; A=[],b=[],Aeq = [],beq = [],lb = [],ub = [],
  usesolver = IpoptSolver(), starter = zeros(size(H,2)))

  r,s = size(H)

  m = Model(solver = usesolver)

  #Initialize the point, you can specify the initial point as you want
  #using the starter input on the function.

  if starter == zeros(size(H,2))
    @variable(m,x[1:s])
  else
    @variable(m, x[i=1:s], start = starter[i])
  end

  #Creates inequality constraint, if A or b are different from the empty array

  if !(isempty(A)) || !(isempty(b))
    @constraintref loopConstraint1[1:length(b)]
    for i = 1:length(b)
      loopConstraint1[i] = @constraint(m,dot(A[i,:],x) <= b[i])
    end
  end

  #Creates equality constraint if it exists.

  if !(isempty(Aeq)) || !(isempty(beq))
    @constraintref loopConstraint2[1:length(beq)]
    for i = 1:length(beq)
      loopConstraint2[i] = @constraint(m, dot(Aeq[i,:],x) == beq[i])
    end
  end

  #Creates the lower bound constraint if it exists

  if !(isempty(lb))
    @constraintref loopConstraint3[1:length(lb)]
    for i = 1:length(lb)
      loopConstraint3[i] = @constraint(m,x[i] >= lb[i])
    end
  end

  #Creates the upper bound constraint if it exists

  if !(isempty(ub))
    @constraintref loopConstraint3[1:length(ub)]
    for i = 1:length(ub)
      loopConstraint4[i] = @constraint(m,x[i] <= ub[i])
    end
  end

  #Specify the objective function, just as in Matlab

  @objective(m, Min, 0.5*dot(H*x,x) + dot(f,x))

  #print the problem on the REPL

  print(m)

  #solve the minimization problem

  solve(m)

  #return the best solution given by the solver

  return getvalue(x)
  
end
