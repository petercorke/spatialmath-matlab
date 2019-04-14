% test classdef methods for array

a = class1(1)
thing(a)

b(1) = class1(2);
b(2) = class1(3);
b(3) = class1(4);
b
fprintf('--- thing(b)\n')
thing(b)
fprintf('--- b.thing\n')
b.thing
fprintf('--- b.thing()\n')
b.thing()
fprintf('---\n')
