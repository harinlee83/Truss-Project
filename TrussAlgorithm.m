[jointNum, memberNum] = size(C);
Cx = zeros(jointNum, memberNum);
Cy = zeros(jointNum, memberNum);
[jointLocation, memberLocation] = find(C == 1);
 
for i = 1:(length(jointLocation)/2)
   vec = find(C(:,i) == 1);
   Cx(jointLocation(2*i-1), memberLocation(2*i-1)) = (X(vec(2))-X((vec(1))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
   Cx(jointLocation(2*i), memberLocation(2*i)) = (X(vec(1)) -X((vec(2))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
end
for i = 1:(length(jointLocation)/2)
   vec = find(C(:,i) == 1);
   Cy(jointLocation(2*i-1), memberLocation(2*i-1)) = (Y(vec(2))-Y((vec(1))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
   Cy(jointLocation(2*i), memberLocation(2*i)) = (Y(vec(1)) -Y((vec(2))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2)); 
end
A = [Cx Sx;Cy Sy];
 
T = A\L;
fprintf('EK301, Section A3, Group: Harin Lee, Landon Kushimi, Brian Mahabir, 4/3/2020 \n');
fprintf('Load: %.6f\n',L(find(L ~= 0)));
fprintf('Member forces in Newtons:\n');
for i = 1:memberNum
    if round(T(i),3) > 0
        fprintf('m%d: %.3f (T)\n',i,abs(T(i)));
    elseif round(T(i),3) == 0
        fprintf('m%d: Zero force member\n',i);
    else
        fprintf('m%d: %.3f (C)\n',i,abs(T(i)));
    end
end
fprintf('Reaction forces in Newtons:\n');
fprintf('Sx1: %.3f \n',round(abs(T(end-2)),2,'significant'));
fprintf('Sy1: %.3f \n',round(T(end-1),2,'significant'));
fprintf('Sy2: %.3f \n',round(T(end),2,'significant'));
totalLength = 0;
for i = 1:memberNum
            vec = find(C(:,i) == 1);
           %Uncomment to see the lengths fprintf('Member %d length: %f cm\n',i,sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2)));
            totalLength = totalLength + sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
end
cost = 10 * jointNum + totalLength;
fprintf('Cost of truss: $%.2f\n', cost);
 
%%Uncertainty code using Euler fit
Buckling_load = zeros(memberNum,1);
for i = 1:memberNum
            vec = find(C(:,i) == 1);
            memberlength = sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
            %%fprintf('Member length %d is %f \n',i,memberlength)%%
            Upper_Confidence = 1400/(memberlength.^2) + 2*1.6; 
            Lower_Confidence = 1400/(memberlength.^2) - 2*1.6;
            Buckling_load(i) = (Upper_Confidence+Lower_Confidence)/2;
            uncertainty_member_percentage(i) = (Buckling_load(i)-Lower_Confidence)./Buckling_load(i); %%find uncertainty as a
            %percentage
end

%{
for i=1:memberNum
    if round(T(i),3) <  0 
      fprintf('Buckling load for member %d is %.3f \n',i,Buckling_load(i));
    else
        fprintf('N/A\n');
    end
end
%}

%%resize T into a new matrix to not include support forces
newL = L;
newL(find(L ~= 0)) = 1; %Sets load to 1N in order to find SR
oneNewtonT = A\newL; %T matrix where load = 1N
newT = oneNewtonT(1:length(oneNewtonT)-3,1); %%delete last 3 forces which are support 
 
%%find SR using the buckling force
SR = length(newT); %allocate SR with all members in vector including the tension members and zero force members
 
for i=1:length(newT)
    if newT(i) < 0 
      SR(i) = abs(newT(i)./Buckling_load(i));
    else
        SR(i) = 0; %For zero force members and tension members, set SR to zero
    end
end
 
%Calculate failure
[SR_max, idx] = max(SR(:)); %make sure to get index of max SR
fprintf('Member %d is the critical member (and all others like it)\n',idx)
maxFailureLoad = 1/SR_max;
max_uncertainty = uncertainty_member_percentage(idx)*maxFailureLoad; %use index to find 
% fprintf('Member %d will fail first\n',find(SR == max(SR(:))));
fprintf('The max theoretical load is: %f N\n',maxFailureLoad);
fprintf('The uncertainty in that max theoretical load is : %f N\n',max_uncertainty);
maxRatio = maxFailureLoad/cost;
fprintf('Theoretical max load/cost ratio in N/$: %.5f\n', maxRatio);
