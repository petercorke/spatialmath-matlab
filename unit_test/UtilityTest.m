%% This is for testing the Utility functions in the robotics Toolbox
function tests = UtilityTest
  tests = functiontests(localfunctions);
end


%    circle                     - compute/draw points on a circle
function circle_test(tc)
    C = [1 2]';
    R = 5;
    x = circle(C, R, 'n', 15 );
    tc.verifyClass(x, 'double');
    tc.verifySize(x, [2 15]);
    tc.verifyTrue( all( (colnorm(x - C(:)) - R) < 1e-10) )
    
    circle(C, R, 'n', 15 );
    
    C = [1 2];
    x = circle(C, R, 'n', 15 );
    tc.verifyClass(x, 'double');
    tc.verifySize(x, [2 15]);
    tc.verifyTrue( all( (colnorm(x - C(:)) - R) < 1e-10) )
    
    circle(C, R, 'n', 15 );
    
    C = [1 2 3];
    x = circle(C, R, 'n', 15 );
    tc.verifyClass(x, 'double');
    tc.verifySize(x, [3 15]);
    tc.verifyTrue( all( (colnorm(x(1:2,:) - C(1:2)') - R) < 1e-10) )
    tc.verifyTrue( all( abs(x(3,:) - C(3)) < 1e-10) )
    
    circle(C, R, 'n', 15 );
   
end

%    colnorm                    - columnwise norm of matrix
function co_testlnorm(tc)
    x= [6.0000    2.5451   -3.0451   -3.0451    2.5451
        2.0000    6.7553    4.9389   -0.9389   -2.7553];
    cn = colnorm(x);
    tc.verifyEqual(cn, [6.3246    7.2188    5.8022    3.1866    3.7509],...
                                  'absTol',1e-4);
end
    
              
%    isvec                      - true if argument is a 3-vector
function is_testvec(tc)
    vh = [1 2 3];
    vv = [1;2;3];
    s = 45;
    tc.verifyTrue(isvec(vh),3);
    tc.verifyTrue(isvec(vv));
    verifyFalse(tc, isvec(s));
    verifyFalse(tc, isvec(ones(2,2)));
    verifyFalse(tc, isvec(ones(2,2,2)));
end
    
%    numcols                    - number of columns in matrix
function numcols_test(tc)
    a = ones(2,3,4);
    b = 2;
    tc.verifyEqual(numcols(a),3);
    tc.verifyEqual(numcols(b),1);
end

%    numrows                    - number of rows in matrix
function numrows_test(tc)
    a = ones(2,3,4);
    b = 3;
    tc.verifyEqual(numrows(a),2);
    tc.verifyEqual(numrows(b),1);
end

%    Polygon                    - general purpose polygon class
function Po_testlygon(tc)
    v = [1 2 1 2;1 1 2 2];
    p = Polygon(v);
%    unit                       - unitize a vector
end

function unit_test(tc)
    vh = [1 2 3];
    vv = [1;2;3];
    vo = [0 0 0];
    tc.verifyEqual(unit(vh), [0.2673    0.5345    0.8018], 'absTol',1e-4);
    tc.verifyEqual(unit(vv), [0.2673
                                         0.5345
                                         0.8018], 'absTol',1e-4);

    verifyError(tc,  @() unit(vo), 'SMTB:unit:zero_norm');
end

%    tb_optparse                - toolbox argument parser
function tb_optparse_test(tc)

    opt.one = [];
    opt.two = 2;
    opt.three = 'three';
    opt.four = false;
    opt.five = true;
    opt.color = {'red', 'green', 'blue'};
    opt.select = {'#bob', '#nancy'};

    opt2 = tb_optparse(opt, {'verbose'});
    tc.verifyEqual(opt2.one, []);
    tc.verifyEqual(opt2.two, 2);
    tc.verifyEqual(opt2.three, 'three');
    tc.verifyEqual(opt2.four, false);
    tc.verifyEqual(opt2.five, true);
    tc.verifyEqual(opt2.color, 'red');
    tc.verifyEqual(opt2.select, 1);

    opt2 = tb_optparse(opt, {'one', 7});
    tc.verifyEqual(opt2.one, 7);
    tc.verifyEqual(opt2.two, 2);
    tc.verifyEqual(opt2.three, 'three');
    tc.verifyEqual(opt2.four, false);
    tc.verifyEqual(opt2.five, true);
    tc.verifyEqual(opt2.color, 'red');
    tc.verifyEqual(opt2.select, 1);

    verifyError(tc, @() tb_optparse(opt, {'one'}), 'SMTB:tboptparse:badargs');

    opt2 = tb_optparse(opt, {'two', 3});
    tc.verifyEqual(opt2.one, []);
    tc.verifyEqual(opt2.two, 3);
    tc.verifyEqual(opt2.three, 'three');
    tc.verifyEqual(opt2.four, false);
    tc.verifyEqual(opt2.five, true);
    tc.verifyEqual(opt2.color, 'red');
    tc.verifyEqual(opt2.select, 1);

    opt2 = tb_optparse(opt, {'three', 'bob'});
    tc.verifyEqual(opt2.one, []);
    tc.verifyEqual(opt2.two, 2);
    tc.verifyEqual(opt2.three, 'bob');
    tc.verifyEqual(opt2.four, false);
    tc.verifyEqual(opt2.five, true);
    tc.verifyEqual(opt2.color, 'red');
    tc.verifyEqual(opt2.select, 1);

    opt2 = tb_optparse(opt, {'four'});
    tc.verifyEqual(opt2.one, []);
    tc.verifyEqual(opt2.two, 2);
    tc.verifyEqual(opt2.three, 'three');
    tc.verifyEqual(opt2.four, true);
    tc.verifyEqual(opt2.five, true);
    tc.verifyEqual(opt2.color, 'red');
    tc.verifyEqual(opt2.select, 1);

    opt2 = tb_optparse(opt, {'nofour'});
    tc.verifyEqual(opt2.one, []);
    tc.verifyEqual(opt2.two, 2);
    tc.verifyEqual(opt2.three, 'three');
    tc.verifyEqual(opt2.four, false);
    tc.verifyEqual(opt2.five, true);
    tc.verifyEqual(opt2.color, 'red');
    tc.verifyEqual(opt2.select, 1);

    opt2 = tb_optparse(opt, {'nofive'});
    tc.verifyEqual(opt2.one, []);
    tc.verifyEqual(opt2.two, 2);
    tc.verifyEqual(opt2.three, 'three');
    tc.verifyEqual(opt2.four, false);
    tc.verifyEqual(opt2.five, false);
    tc.verifyEqual(opt2.color, 'red');
    tc.verifyEqual(opt2.select, 1);

    opt2 = tb_optparse(opt, {'green'});
    tc.verifyEqual(opt2.one, []);
    tc.verifyEqual(opt2.two, 2);
    tc.verifyEqual(opt2.three, 'three');
    tc.verifyEqual(opt2.four, false);
    tc.verifyEqual(opt2.five, true);
    tc.verifyEqual(opt2.color, 'green');
    tc.verifyEqual(opt2.select, 1);

    opt2 = tb_optparse(opt, {'nancy'});
    tc.verifyEqual(opt2.one, []);
    tc.verifyEqual(opt2.two, 2);
    tc.verifyEqual(opt2.three, 'three');
    tc.verifyEqual(opt2.four, false);
    tc.verifyEqual(opt2.five, true);
    tc.verifyEqual(opt2.color, 'red');
    tc.verifyEqual(opt2.select, 2);

    opt2 = tb_optparse(opt, {});
    tc.verifyEqual(opt2.verbose, false);
    tc.verifyEqual(opt2.debug, 0);
    opt2 = tb_optparse(opt, {'verbose'});
    tc.verifyEqual(opt2.verbose, true);
    opt2 = tb_optparse(opt, {'verbose=2'});
    tc.verifyEqual(opt2.verbose, 2);
    opt2 = tb_optparse(opt, {'debug', 11});
    tc.verifyEqual(opt2.debug, 11);

    opt2 = tb_optparse(opt, {'showopt'});

    opt3.color = 'green';
    opt3.five = false;

    opt2 = tb_optparse(opt, {'setopt', opt3});
    tc.verifyEqual(opt2.one, []);
    tc.verifyEqual(opt2.two, 2);
    tc.verifyEqual(opt2.three, 'three');
    tc.verifyEqual(opt2.four, false);
    tc.verifyEqual(opt2.five, false);
    tc.verifyEqual(opt2.color, 'green');
    tc.verifyEqual(opt2.select, 1);

    [opt2,args] = tb_optparse(opt, {1, 'three', 4, 'spam', 2, 'red',  'spam'});
    tc.verifyEqual(args, {1, 'spam', 2, 'spam'});

    verifyError(tc,  @() tb_optparse(opt, {'two'}), 'SMTB:tboptparse:badargs');
    verifyError(tc,  @() tb_optparse(opt, 'bob'), 'SMTB:tboptparse:badargs');
end


function Animate_folder_test(tc)
    file = 'test';
    a = Animate(file);
    plot([1 2], [3 4]);
    for i=1:10
        a.add();
    end
    a.close();
    
    tc.verifyTrue( exist(file, 'dir') > 0 );
    
    d = dir('test/*.png');
    tc.verifyEqual( length(d), 10);
    
    for dd=d'
        delete(fullfile(file, dd.name));
    end
    rmdir(file)
end

function Animate_mp4_test(tc)
    tc.assumeTrue(ismac || ispc);  % this test won't work on linux
    
    file = 'test.mp4';
    a = Animate(file);
    plot([1 2], [3 4]);
    for i=1:10
        a.add();
    end
    a.close();
    
    tc.verifyTrue( exist(file, 'file') > 0 );
    
    delete(file)
end

function Animate_gif_test(tc)
    file = 'test.gif';
    a = Animate(file);
    plot([1 2], [3 4]);
    for i=1:10
        a.add();
    end
    a.close();
    
    tc.verifyTrue( exist(file, 'file') > 0);
    
    im = imread(file);
    tc.verifyEqual( size(im,4), 10);
        delete(file)
end

function homline_test(tc)
    L = homline(-2, 3, 3, 3);  % y = 3
    tc.verifyEqual( L/L(2), [0 1 -3]); 
    
    L = homline(-2, 10, -2, 30);   % x = -2
    tc.verifyEqual( L/L(1), [1 0 2]);
    
    L = homline(1, 5, 3, 8);   % x = -2
    tc.verifyEqual( L/L(1)*3, [3 -2 7]);
end

function about_test(tc)
    
    a  = 3;
    
    about(a)
    about('a')  % command line equivalent
    
    tc.verifyError( @() about('b'), 'SMTB:about')
    
    a = sqrt(-1)
    about(a)
    
    a = rand(10,10);
    about(a)
    
    a = SE3();
    about(a)

    a = zeros(1000000,1);
    about(a)
end

function randinit_test(tc)
    
    randinit()
    x = rand;
    rand(100,1);
    randinit()
    y = rand;
    
    tc.verifyEqual(x, y);
end

function angdiff_test(tc)
    % 2 arg case
    tc.verifyEqual( angdiff(2, 1), 1);
    tc.verifyEqual( angdiff(1, 2), -1);
    tc.verifyEqual( angdiff(2, 2), 0);
    
    % 1 arg case
    tc.verifyEqual( angdiff([2, 1]), 1);
    tc.verifyEqual( angdiff([1, 2]), -1);
    tc.verifyEqual( angdiff([2, 2]), 0);
    
    pi34 = pi*3/4;
    pi2 = pi/2;
    % 2 arg case
    tc.verifyEqual( angdiff(pi34, pi34), 0)
    tc.verifyEqual( angdiff(pi34, -pi34), -pi2)
    tc.verifyEqual( angdiff(-pi34, pi34), pi2)
    
    % 1 arg case
    tc.verifyEqual( angdiff([pi34, pi34]), 0)
    tc.verifyEqual( angdiff([pi34, -pi34]), -pi2)
    tc.verifyEqual( angdiff([-pi34, pi34]), pi2)
    
    % vector case
    tc.verifyEqual( angdiff( [pi34 pi34 -pi34], [pi34 -pi34 pi34]), [0 -pi2 pi2])
    
    tc.verifyError( @() angdiff( [pi34 pi34 -pi34], [pi34 -pi34]), 'SMTB:angdiff:badarg')
end

function stl_test(tc)
    [v, f, n, name] = stlRead('data/20mm_cube_ascii.stl');
    tc.verifyClass(v, 'double');
    tc.verifySize(v, [8 3]);
    tc.verifyClass(f, 'double');
    tc.verifySize(f, [12 3]);
    tc.verifyClass(n, 'double');
    tc.verifySize(n, [12 3]);
    
    [v, f, n, name] = stlRead('data/20mm_cube_binary.stl');
    tc.verifyClass(v, 'double');
    tc.verifySize(v, [8 3]);
    tc.verifyClass(f, 'double');
    tc.verifySize(f, [12 3]);
    tc.verifyClass(n, 'double');
    tc.verifySize(n, [12 3]);
end
