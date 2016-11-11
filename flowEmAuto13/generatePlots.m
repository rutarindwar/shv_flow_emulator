#!/usr/bin/octave -qf


load ("anggrid.txt");
load("vgrid.txt");
load("ypredb.txt");
load("vcontoursaveb.txt");
x = load("designp_000250.txt");
y = load("logheightp_000250.txt");
size(vcontoursaveb)
%figure(1)
%s1 = surf(vgrid, anggrid, ypredb');
%shading flat
%set(s1,'facealpha',0.2)
%hold on
%plot3(x(:, 1), x(:, 2)*pi/180, y, 'k*');
plot(anggrid,vcontoursaveb);
print("Contourlines.pdf", "-dpdf");
pause()
