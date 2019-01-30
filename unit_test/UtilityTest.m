%% This is for testing the Utility functions in the robotics Toolbox
function tests = UtilityTest
  tests = functiontests(localfunctions);
end

%% Utility
%    about                      - summary of object size and type
function ab_testout(testCase)
    a = [1 2 3; 4 5 6];
    about(a);
    about a;
end

%    angdiff                    - subtract 2 angles modulo 2pi
function an_testgdiff(testCase)
    verifyEqual(testCase, angdiff(pi/2,pi), -1.5708, 'absTol',1e-4);
end

%    circle                     - compute/draw points on a circle
function ci_testrcle(testCase)
    x= circle([1 2],5,'n', 5 );
    verifyEqual(testCase, x, [6.0000    2.5451   -3.0451   -3.0451    2.5451
                                  2.0000    6.7553    4.9389   -0.9389   -2.7553],...
                                  'absTol',1e-4);
end

%    colnorm                    - columnwise norm of matrix
function co_testlnorm(testCase)
    x= [6.0000    2.5451   -3.0451   -3.0451    2.5451
        2.0000    6.7553    4.9389   -0.9389   -2.7553];
    cn = colnorm(x);
    verifyEqual(testCase, cn, [6.3246    7.2188    5.8022    3.1866    3.7509],...
                                  'absTol',1e-4);
end
    
              
%    isvec                      - true if argument is a 3-vector
function is_testvec(testCase)
    vh = [1 2 3];
    vv = [1;2;3];
    s = 45;
    verifyTrue(testCase, isvec(vh),3);
    verifyTrue(testCase, isvec(vv));
    verifyFalse(testCase, isvec(s));
    verifyFalse(testCase, isvec(ones(2,2)));
    verifyFalse(testCase, isvec(ones(2,2,2)));
end
    
%    numcols                    - number of columns in matrix
function numcols_test(testCase)
    a = ones(2,3,4);
    b = 2;
    verifyEqual(testCase, numcols(a),3);
    verifyEqual(testCase, numcols(b),1);
end

%    numrows                    - number of rows in matrix
function numrows_test(testCase)
    a = ones(2,3,4);
    b = 3;
    verifyEqual(testCase, numrows(a),2);
    verifyEqual(testCase, numrows(b),1);
end

%    Polygon                    - general purpose polygon class
function Po_testlygon(testCase)
    v = [1 2 1 2;1 1 2 2];
    p = Polygon(v);
%    unit                       - unitize a vector
end

function unit_test(testCase)
    vh = [1 2 3];
    vv = [1;2;3];
    vo = [0 0 0];
    verifyEqual(testCase, unit(vh), [0.2673    0.5345    0.8018], 'absTol',1e-4);
    verifyEqual(testCase, unit(vv), [0.2673
                                         0.5345
                                         0.8018], 'absTol',1e-4);

    verifyError(testCase,  @() unit(vo), 'RTB:unit:zero_norm');
end

%    tb_optparse                - toolbox argument parser
function tb_optparse_test(testCase)

    opt.one = [];
    opt.two = 2;
    opt.three = 'three';
    opt.four = false;
    opt.five = true;
    opt.color = {'red', 'green', 'blue'};
    opt.select = {'#bob', '#nancy'};

    opt2 = tb_optparse(opt, {'verbose'});
    verifyEqual(testCase, opt2.one, []);
    verifyEqual(testCase, opt2.two, 2);
    verifyEqual(testCase, opt2.three, 'three');
    verifyEqual(testCase, opt2.four, false);
    verifyEqual(testCase, opt2.five, true);
    verifyEqual(testCase, opt2.color, 'red');
    verifyEqual(testCase, opt2.select, 1);

    opt2 = tb_optparse(opt, {'one', 7});
    verifyEqual(testCase, opt2.one, 7);
    verifyEqual(testCase, opt2.two, 2);
    verifyEqual(testCase, opt2.three, 'three');
    verifyEqual(testCase, opt2.four, false);
    verifyEqual(testCase, opt2.five, true);
    verifyEqual(testCase, opt2.color, 'red');
    verifyEqual(testCase, opt2.select, 1);

    verifyError(testCase, @() tb_optparse(opt, {'one'}), 'RTB:tboptparse:badargs');

    opt2 = tb_optparse(opt, {'two', 3});
    verifyEqual(testCase, opt2.one, []);
    verifyEqual(testCase, opt2.two, 3);
    verifyEqual(testCase, opt2.three, 'three');
    verifyEqual(testCase, opt2.four, false);
    verifyEqual(testCase, opt2.five, true);
    verifyEqual(testCase, opt2.color, 'red');
    verifyEqual(testCase, opt2.select, 1);

    opt2 = tb_optparse(opt, {'three', 'bob'});
    verifyEqual(testCase, opt2.one, []);
    verifyEqual(testCase, opt2.two, 2);
    verifyEqual(testCase, opt2.three, 'bob');
    verifyEqual(testCase, opt2.four, false);
    verifyEqual(testCase, opt2.five, true);
    verifyEqual(testCase, opt2.color, 'red');
    verifyEqual(testCase, opt2.select, 1);

    opt2 = tb_optparse(opt, {'four'});
    verifyEqual(testCase, opt2.one, []);
    verifyEqual(testCase, opt2.two, 2);
    verifyEqual(testCase, opt2.three, 'three');
    verifyEqual(testCase, opt2.four, true);
    verifyEqual(testCase, opt2.five, true);
    verifyEqual(testCase, opt2.color, 'red');
    verifyEqual(testCase, opt2.select, 1);

    opt2 = tb_optparse(opt, {'nofour'});
    verifyEqual(testCase, opt2.one, []);
    verifyEqual(testCase, opt2.two, 2);
    verifyEqual(testCase, opt2.three, 'three');
    verifyEqual(testCase, opt2.four, false);
    verifyEqual(testCase, opt2.five, true);
    verifyEqual(testCase, opt2.color, 'red');
    verifyEqual(testCase, opt2.select, 1);

    opt2 = tb_optparse(opt, {'nofive'});
    verifyEqual(testCase, opt2.one, []);
    verifyEqual(testCase, opt2.two, 2);
    verifyEqual(testCase, opt2.three, 'three');
    verifyEqual(testCase, opt2.four, false);
    verifyEqual(testCase, opt2.five, false);
    verifyEqual(testCase, opt2.color, 'red');
    verifyEqual(testCase, opt2.select, 1);

    opt2 = tb_optparse(opt, {'green'});
    verifyEqual(testCase, opt2.one, []);
    verifyEqual(testCase, opt2.two, 2);
    verifyEqual(testCase, opt2.three, 'three');
    verifyEqual(testCase, opt2.four, false);
    verifyEqual(testCase, opt2.five, true);
    verifyEqual(testCase, opt2.color, 'green');
    verifyEqual(testCase, opt2.select, 1);

    opt2 = tb_optparse(opt, {'nancy'});
    verifyEqual(testCase, opt2.one, []);
    verifyEqual(testCase, opt2.two, 2);
    verifyEqual(testCase, opt2.three, 'three');
    verifyEqual(testCase, opt2.four, false);
    verifyEqual(testCase, opt2.five, true);
    verifyEqual(testCase, opt2.color, 'red');
    verifyEqual(testCase, opt2.select, 2);

    opt2 = tb_optparse(opt, {});
    verifyEqual(testCase, opt2.verbose, false);
    verifyEqual(testCase, opt2.debug, 0);
    opt2 = tb_optparse(opt, {'verbose'});
    verifyEqual(testCase, opt2.verbose, true);
    opt2 = tb_optparse(opt, {'verbose=2'});
    verifyEqual(testCase, opt2.verbose, 2);
    opt2 = tb_optparse(opt, {'debug', 11});
    verifyEqual(testCase, opt2.debug, 11);

    opt2 = tb_optparse(opt, {'showopt'});

    opt3.color = 'green';
    opt3.five = false;

    opt2 = tb_optparse(opt, {'setopt', opt3});
    verifyEqual(testCase, opt2.one, []);
    verifyEqual(testCase, opt2.two, 2);
    verifyEqual(testCase, opt2.three, 'three');
    verifyEqual(testCase, opt2.four, false);
    verifyEqual(testCase, opt2.five, false);
    verifyEqual(testCase, opt2.color, 'green');
    verifyEqual(testCase, opt2.select, 1);

    [opt2,args] = tb_optparse(opt, {1, 'three', 4, 'spam', 2, 'red',  'spam'});
    verifyEqual(testCase, args, {1, 'spam', 2, 'spam'});

    verifyError(testCase,  @() tb_optparse(opt, {'two'}), 'RTB:tboptparse:badargs');
    verifyError(testCase,  @() tb_optparse(opt, 'bob'), 'RTB:tboptparse:badargs');
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
    
    d = dir('test/*.png')
    tc.verifyEqual( length(d), 10);
    
    for dd=d'
        delete(fullfile(file, dd.name));
    end
    rmdir(file)
end

% function Animate_mp4_test(tc)
%     file = 'test.mp4';
%     a = Animate(file);
%     plot([1 2], [3 4]);
%     for i=1:10
%         a.add();
%     end
%     a.close();
%     
%     tc.verifyTrue( exist(file, 'file') > 0 );
%     
%     delete(file)
% end

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
    
    tc.verifyError( @() about('b'), 'RTB:about')
    
    a = sqrt(-1)
    about(a)
    
    a = rand(10,10);
    about(a)
    
    a = SE3();
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
    tc.verifyEqual( angdiff(2, 1), 1);
    tc.verifyEqual( angdiff(1, 2), -1);
    tc.verifyEqual( angdiff(2, 2), 0);
    
    pi34 = pi*3/4;
    pi2 = pi/2;
    tc.verifyEqual( angdiff(pi34, pi34), 0)
    tc.verifyEqual( angdiff(pi34, -pi34), -pi2)
    tc.verifyEqual( angdiff(-pi34, pi34), pi2)
end