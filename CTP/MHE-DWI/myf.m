function [ fg ] = myf( gg, ug )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global Dact t Sd Wm

util = T_G_Utility(Sd(:, t), Dact(t), gg, Wm, t);

fg = util - ug;

end

