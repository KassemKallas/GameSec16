function [v,p,q] = game_solve(A)
% Solve a zero-sum matrix game A for the row player.
%
% The input to this function is a matrix A representing a classical matrix
% game in which A(i,j) is the payoff from the column player to the row
% player. 
%
% The outputs are [v,p] where
% * v is the value of the game and
% * p is the mixed strategy for the row player.
% 
% To calculate the mixed strategy for the column player, call
% game_solve(-A'). Alternatively, call
% [v,p,q] = game_solve(A)
% and q will be the mixed strategy for the column player.
%
% Example: The matrix for rock-paper-scissors is
% A = [0 -1 1; 1 0 -1; -1 1 0]
%
% calling [v,p] = game_solve(A) yields
%
% v =
%    0
% p =
%    0.3333
%    0.3333
%    0.3333
%
% showing that one should play rock, paper, or scissors each with
% probability 1/3 and the value of the game is zero.

[r,c] = size(A);

AA = [-A', ones(c,1)];
Aeq = [ones(1,r),0];

b = zeros(c,1);
beq = 1;
lb = [zeros(r,1);-inf];
f = [ zeros(r,1);-1];
options = optimset('Display', 'off');

p = linprog(f,AA,b,Aeq,beq,lb,[],[],options);
v = p(r+1);
p = p(1:r);

if nargout > 2
    [w,q] = game_solve(-A');
end
