% Translate data from student input format to FEM format required 
% for FEM solution.  
% 	Paul E. Barbone			April 2020 
%
% Declare global variables:=======================================
%
%	Nel:		Total number of elements
%	NNodes: 	Total number of nodes.
%	ElType: 	Type of elements used to make up structure.  
%			 = 1: all (2D) bar elements in uniaxial loading.
%			 = 2: all (2D) beam elements with both bending 
%				and elongation.  
%			 = 3: all (2D) bar elements with buckling. 
%	XNodes:		Nnodes x 2 matrix containing (x,y) coordinates 
%			for each node.  This information must be supplied 
%			by the user in "input.m"
%	ElNodes: 	Nel x 2 matrix containing the node numbers for 
%			each element.  
%	LtoG:		Nel x (2Ndof) matrix containing the global eqn
%			number for each local eqn number.  This matrix is 
%			automatically generated in function "MakeLtoG" 
%			according to the convention described there.
%	BigK		Global stiffness matrix.  Boundary values have 
%			not been taken into account in forming BigK.  
%	BigF		The global matrix of forces applied to the structure.
%			See the description of the equation ordering in 
%			function "MakeLtoG", and notes below.
% 	BCs 		(NBCs x 2) array, where NBCs = number of known 
% 			displacement components. The first column in BCs 
% 			is P, the eqn number for the known dof.  The 
% 			second column contains the known value of the dof.  
%	ElProps		Nel x ...  Physical properties of each element. 
%			The size of the matrix (the number of columns) 
%			depends on the number of physical properties required
%			to characterize that type of element.  
%
% LOCAL VARIABLES (in this function only) 
%
%	Ndof: 		Number of dof (degrees of freedom or 
%			unknowns) per node.

global nel nnodes eltype bigXnodes xnodes elnodes ltog bigk bigf bcs elprops;

% ========================================================================
% BEGIN USER INPUT:
% From the student input file, we have: 
					% C = Connectivity matrix
					% L = a load vector, determines loaded node.
					% X = nodal coordinate X
					% Y = nodal coordinate Y
					% Sx = determines fixed node
					% Sy = determines roller node
[nnodes, nel] = size(C); 
 		% nel = Total number of elements in problem. 
		% nnodes = Total number of nodes.
% END USER INPUT
% ========================================================================

% ========================================================================
% Make XNodes array: 
% Check that X and Y have same length: 
if (length(X) - nnodes) ~=0 
  disp('X and C have mismatched dimensions') 
  return
end
if (length(Y) - nnodes) ~=0 
  disp('Y and C have mismatched dimensions') 
  return
end
if (length(Y) - length(X)) ~=0 
  disp('X and Y have mismatched dimensions') 
  return
end
% Create bigXnodes: 
bigXnodes(:,1) = X;  
bigXnodes(:,2) = Y;  

% Check dimensions of truss: 
if (max(Y) > 0)
  disp('Truss should lie below the line y= 0 .') 
  return
end
if (abs(min(X)) > 1e-12 ) 
  disp('Left end of truss should be at location x=0.' ) 
  return
end
if (abs(max(X)-56) > 1e-12 ) 
  disp(['Right end of truss should be at location x=56cm.  Are distances entered in cm?' ]) 
  return
end


% ========================================================================
% Make ElNodes array:
elnodes = zeros(nel, 2);	% Initialize ElNodes for elements with 2 nodes
% Check connectivity matrix so that each element has only two nodes.  
check = sum(C,1); 	% Sum each column of C to get number of nodes/element.
if (check ~= 2*ones(size(check))) 
   disp(['C shows some members connected to more than two joints.']) 
   return
end

for jel = 1:nel
  tn = find (C(:,jel) ==1) ; 	% find nodes connected to element jel
  elnodes(jel,:) = [tn(1),tn(2)]; 		% put them in elnodes.
end

% ========================================================================
% Compute length of each element/member and check that it is within range: 
for n = 1:nel
   x1 = bigXnodes(elnodes(n,1),:);
   x2 = bigXnodes(elnodes(n,2),:);
   bigL(n) = norm(x2-x1); 
   if or((bigL(n) < 9),(bigL(n)>16))
     disp(['The length of the member connecting joints ', num2str(elnodes(n,1)), ' and ', num2str(elnodes(n,2)), ' is ', num2str(bigL(n)), 'cm, which is outside of range.'])
     % return
   end
end

%
% ========================================================================
% Define Loading.  "Loading" means two things:  prescribed forces and 
% prescribed boundary conditions.   We need to enter both.  
%
% Find loaded node from L = (2*NNodes) x 1: 
% The numbering used in L is:
% the first NNodes are eqns in the x-direction
% the second set of NNodes are eqns in the y-direction
loaded_node = find(L~= 0 ); 
if (length(loaded_node) ~=1)
  disp('Something is wrong with L.') 
  return
end
loaded_node = loaded_node - nnodes; 	% subtract x equations to get 
					% which joint is loaded in y-dir. 
if (loaded_node < 0 )
  disp('Something is wrong with L.') 
  return
elseif (loaded_node > nnodes )
  disp('Something is wrong with L.') 
  return
end
% Check that loaded node is within the prescribed tolerance: 
if or((bigXnodes(loaded_node,1) > 27), (bigXnodes(loaded_node,1) < 25))
  disp(['The loaded joint ,', num2str(loaded_node),', is at location ', num2str(bigXnodes(loaded_node,1)), ' cm.']) 
  disp('Loaded joint location is not within 25-27cm') 
  % return
end 



% ========================================================================
% Determine those nodes that are fixed in each direction: 

% Check Sx: 
check = sum(sum(Sx));
if check ~=1 
  disp('Unexpected structure of Sx');
  return
end
cols = sum(Sx,1);  
col = find(cols ==1);
xy_fixed_node = find(Sx(:,col)==1); 

% Check Sy: 
if Sy(xy_fixed_node,1 +mod(col,3)) ~=1 
  disp('Unexpected structure of Sy');
  return
end
y_fixed_node = find(Sy(:,1+mod(col+1,3))==1); 

trusslength = abs(bigXnodes(y_fixed_node,1) - bigXnodes(xy_fixed_node,1)); 
if (abs(trusslength - 56) > 1e-12) 
  disp('distance between supports is not 56 cm?');  
  disp('Check Sx and Sy');  
end

% ========================================================================
% Draw truss schematic for inspection by student: 
drawschematic(bigXnodes, elnodes,loaded_node, xy_fixed_node, y_fixed_node);     







