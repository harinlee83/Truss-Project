  function drawschematic(xnodes, elnodes, fnode, pinned_node, roller_node);     

% Variables:=======================================
%
%	Nel:		Total number of elements
%	NNodes: 	Total number of nodes.
%	XNodes:		Nnodes x 2 matrix containing (x,y) coordinates 
%			for each node.  This information must be supplied 
%			by the user in "input.m"
%	ElNodes: 	Nel x 2 matrix containing the node numbers for 
%			each element.  
% ========================================================================


nel = size(elnodes,1);	 % number of elements
nnodes = size(xnodes,1); % number of nodes.
			

H = figure(314); 	% Make figure. 

plot(xnodes(:,1), xnodes(:,2), 'ro');
			% plot nodes with red circles.
hold on; 		% keep the nodes in place while we add the elements.
% Add triangle to pinned node: 
plot( xnodes(pinned_node,1), xnodes(pinned_node,2)+1, 'kv', 'MarkerSize',14, 'LineWidth', 3)
% Add large circle to node with roller support: 
plot( xnodes(roller_node,1), xnodes(roller_node,2)+1, 'ko', 'MarkerSize',14, 'LineWidth', 3)

ellabel = char([1:nel] + 64); 	% get character strings for each element number. 
for n = 1:nel		% plot each element one at a time
  n1 = elnodes(n,1);    % Determine the node numbers attached to element.
  n2 = elnodes(n,2);    % Determine the node numbers attached to element.
  x1 = xnodes(n1,1);    % Get (x,y) coordinates of node n1.
  y1 = xnodes(n1,2);    % Get (x,y) coordinates of node n1.
  x2 = xnodes(n2,1);    % Get (x,y) coordinates of node n2.
  y2 = xnodes(n2,2);    % Get (x,y) coordinates of node n2.
  plot([x1,x2],[y1,y2],'b-', 'LineWidth', 2);
			% Draw the element line segment
  text((x1+x2)/2, (y1+y2)/2, ellabel(n), 'FontSize', 14); 
			% Label each element with an element letter.
end;

% Resize plot to give correct proportions to truss: 
xmin = min(xnodes(:,1)); 
xmax = max(xnodes(:,1)); 
ymin = min(xnodes(:,2)); 
ymax = max(xnodes(:,2)); 
brdr = 10; 
axis('equal') 
% axis([ xmin - hght/2, xmax + hght/2, ymin - hght/2, ymax + hght/2]) 
axis([ xmin - brdr, xmax + brdr, ymin - brdr, ymax + brdr]) 

% Label each node with its node number: 
for inode = 1:nnodes
  text(xnodes(inode,1)+1, xnodes(inode,2), num2str(inode), 'FontSize', 14)
end

% Add an arrow to the node with the force: 
quiver( xnodes(fnode,1), xnodes(fnode,2), 0, -7, 0, 'linewidth',3,'color','k', 'MaxHeadSize', 1.0); 

% Add labels
title({'Truss Schematic'}, 'FontSize', 14)
xlabel('Horizontal Position (cm)', 'FontSize', 14)
ylabel('Vertical Position (cm)', 'FontSize', 14)

return
  






