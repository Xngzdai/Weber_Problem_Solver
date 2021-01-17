function Opt_Location = Weber_Solver(W,X,Y0,err)
%   Weber_Solver Opt_Location = Weber_Solver(W,X,Y0) produces the optimal 
%   solution to the Fermat-Weber porblem in a weighted situation. 
%   This function is a realization of the algorithm written in 
%   'Iterative solutions of the Fermat, Weber and attraction-repulsion 
%   problems', from https://en.wikipedia.org/wiki/Weber_problem.
%   W is a row vector of size n, composed of weights of each point. X is a 
%   2-by-n matrix composed of n 2D column vectors. Y0 is a column vector 
%   of size 2 that serves as an initiating point for the optimization. err
%   is a scalar value determining the precision degree of the optimization.
%   The suggestion of default value is 0.001.
%

    %   Plot the initial situation
    distribution_plotter(W,X);
    %   Initialize the output as an two dimensional column vector.
    Y_next = nan(2,1);
    
    %   Repeated optimization process to calculate the optimal position.
    while true
        Y_next = KK_method(W,X,Y0);
        trace_plotter(Y0,Y_next);
        if loop_condition_check(Y0,Y_next,err)
            break
        end
        Y0 = Y_next;
    end
    
    %   Output final answer.
    format long
    Opt_Location = Y_next;
    
    %   Set the hold state of the graph to off
    hold off
    
end


function distance = get_Distance(X,Y0)
%   GET_DISTANCE distance = get_Distance(X,Y0) calculates the distance
%   between two matrices composed by column vectors. It is different
%   and may not be replaceable with the built-in function dist() because
%   get_Distance processes multiple vectors at one time.
%
    sz = size(X);
    len = sz(2);
    Y0 = repmat(Y0,1,len);
    sq_diff = (X - Y0) .^ 2;
    distance = (sq_diff(1,:) + sq_diff(2,:)).^(0.5);
end

function Y_next = KK_method(W,X,Y0)
%   KK_METHOD [Y_next, Y0] = KK_method(Y_next,Y0) completes the algorithm
%   by computing all ingredients together. This algorithm was proposed by 
%   Kuhn and Kuenne (1962). <Kuhn, Harold W. and Robert E. Kuenne, 
%   1962, "An Efficient Algorithm for the Numerical Solution of the 
%   Generalized Weber Problem in Spatial Economics." 
%   Journal of Regional Science 4, 21–34.>
%
    %   General coefficient
    denominator = get_Distance(X,Y0);
    coeff = 1 / sum(W ./ denominator);
    
     %   Each Term
     x = sum(X(1,:) .* W ./ denominator);
     y = sum(X(2,:) .* W ./ denominator);
     term = [x;y]; 
     
     %   Answer
     Y_next = term .* coeff;
end

function TorF = loop_condition_check(Y0,Y_next,err)
%   LOOP_CONDITION_CHECK TorF = loop_condition_check(Y0,Y_next) 
%   checks if the Y0 has converged properly Y_next. If the difference
%   between Y_next and Y0 can be regarded as infinitesimal (based on the
%   user's standard), Y_next is then good enough to output. Otherwise, this
%   judgement leads to next round optimization.
% 
    Y_diff = get_Distance(Y0,Y_next);
    
    if Y_diff >= err
        TorF = 0;
    else
        TorF = 1;
    end
end

function distribution_plotter(W,X)
%   DISTRIBUTION_PLOTTER distribution_plotter(W,X) maps out the location of
%   initial points on an Euclidean surface. Points are colored according to
%   the weight vector; the closer the color is to the red, the larger
%   weight it has and vice versa.
%
    x = X(1,:);
    y = X(2,:);
    scatter(x,y,20,W,'filled')
    hold on
end

function trace_plotter(p1,p2)
%   TRACE_PLOTTER trace_plotter(p1,p2) traces the direction that
%   optimization is achieved. The vector is colored only for distinguishing
%   one another. Reference: Star Strider (29 Oct 2014)https://www.mathworks
%   .com/matlabcentral/answers/160487-how-can-i-draw-a-line-with-arrow-head
%   -between-2-data-points-in-a-plot#answer_156958.
%
    difference = p2-p1; 
    quiver(p1(1),p1(2),difference(1),difference(2),0)
end