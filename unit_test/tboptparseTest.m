function tests = tboptparseTest()
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
   opt.foo = false;
   opt.bar = true;
   opt.blah = [];
   opt.stuff = {};
   opt.choose = {'this', 'that', 'other'};
   opt.select = {'#no', '#yes'};
   opt.old = '@foo';
   opt.d_3d = false;
   tc.TestData.opt = opt;
end

function boolTest(tc)
    opt.foo = false;
    [out,args] = tb_optparse(opt, {});
    tc.verifyEqual(out.foo, false);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo'});
    tc.verifyEqual(out.foo, true);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'distract'});
    tc.verifyEqual(out.foo, false);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'foo', 'distract'});
    tc.verifyEqual(out.foo, true);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'distract', 'foo'});
    tc.verifyEqual(out.foo, true);
    tc.verifyEqual(args, {'distract'});
    
end

function noboolTest(tc)
    opt.foo = true;
    [out,args] = tb_optparse(opt, {});
    tc.verifyEqual(out.foo, true);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo'});
    tc.verifyEqual(out.foo, true);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'nofoo'});
    tc.verifyEqual(out.foo, false);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'distract'});
    tc.verifyEqual(out.foo, true);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'foo', 'distract'});
    tc.verifyEqual(out.foo, true);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'distract', 'foo'});
    tc.verifyEqual(out.foo, true);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'nofoo', 'distract'});
    tc.verifyEqual(out.foo, false);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'distract', 'nofoo'});
    tc.verifyEqual(out.foo, false);
    tc.verifyEqual(args, {'distract'});
end

function setTest(tc)
    opt.foo = [];
    
    [out,args] = tb_optparse(opt, {});
    tc.verifyEqual(out.foo, []);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo', 3});
    tc.verifyEqual(out.foo, 3);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo', 'bar'});
    tc.verifyEqual(out.foo, 'bar');
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo', [1 2 3]});
    tc.verifyEqual(out.foo, [1 2 3]);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo', {1 2 3}});
    tc.verifyEqual(out.foo, {1 2 3});
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'distract'});
    tc.verifyEqual(out.foo, []);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'distract', 'foo', 'bar'});
    tc.verifyEqual(out.foo, 'bar');
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'foo', 'bar', 'distract'});
    tc.verifyEqual(out.foo, 'bar');
    tc.verifyEqual(args, {'distract'});
end

function cellsetTest(tc)
    opt.foo = {};
    
    [out,args] = tb_optparse(opt, {});
    tc.verifyEqual(out.foo, {});
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo', 3});
    tc.verifyEqual(out.foo, {3});
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo', 'bar'});
    tc.verifyEqual(out.foo, {'bar'});
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo', [1 2 3]});
    tc.verifyEqual(out.foo, {[1 2 3]});
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo', {1 2 3}});
    tc.verifyEqual(out.foo, {1 2 3});
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'distract'});
    tc.verifyEqual(out.foo, {});
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'distract', 'foo', 'bar'});
    tc.verifyEqual(out.foo, {'bar'});
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'foo', 'bar', 'distract'});
    tc.verifyEqual(out.foo, {'bar'});
    tc.verifyEqual(args, {'distract'});
end

function chooseTest(tc)
    opt.choose = {'this', 'that', 'other'};
    
    [out,args] = tb_optparse(opt, {});
    tc.verifyEqual(out.choose, 'this');
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'this'});
    tc.verifyEqual(out.choose, 'this');
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'that'});
    tc.verifyEqual(out.choose, 'that');
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'other'});
    tc.verifyEqual(out.choose, 'other');
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'yetanother'});  % this one not in the list
    tc.verifyEqual(out.choose, 'this'); % return default
    tc.verifyEqual(args, {'yetanother'}); % and the arg is returned here
    
    [out,args] = tb_optparse(opt, {'distract', 'that'});
    tc.verifyEqual(out.choose, 'that');
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'that','distract'});
    tc.verifyEqual(out.choose, 'that');
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'yetanother','distract'});
    tc.verifyEqual(out.choose, 'this');
    tc.verifyEqual(args, {'yetanother','distract'});
end

function hashchooseTest(tc)
    opt.choose = {'#this', '#that', '#other'};
    
    [out,args] = tb_optparse(opt, {});
    tc.verifyEqual(out.choose, 1);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'this'});
    tc.verifyEqual(out.choose, 1);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'that'});
    tc.verifyEqual(out.choose, 2);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'other'});
    tc.verifyEqual(out.choose, 3);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'yetanother'});  % this one not in the list
    tc.verifyEqual(out.choose, 1); % return default
    tc.verifyEqual(args, {'yetanother'}); % and the arg is returned here
    
    [out,args] = tb_optparse(opt, {'distract', 'that'});
    tc.verifyEqual(out.choose, 2);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'that','distract'});
    tc.verifyEqual(out.choose, 2);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'yetanother','distract'});
    tc.verifyEqual(out.choose, 1);
    tc.verifyEqual(args, {'yetanother','distract'});
end

function synonymTest(tc)
    opt.foo = false;
    opt.bar = '@foo';
    
    [out,args] = tb_optparse(opt, {});
    tc.verifyEqual(out.foo, false);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'foo'});
    tc.verifyEqual(out.foo, true);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'bar'});
    tc.verifyEqual(out.foo, true);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'distract'});
    tc.verifyEqual(out.foo, false);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'bar', 'distract'});
    tc.verifyEqual(out.foo, true);
    tc.verifyEqual(args, {'distract'});
    
    [out,args] = tb_optparse(opt, {'distract', 'bar'});
    tc.verifyEqual(out.foo, true);
    tc.verifyEqual(args, {'distract'});
end

function startDigitTest(tc)
    % bool
    opt.d_3foo = false;
    [out,args] = tb_optparse(opt, {});
    tc.verifyEqual(out.d_3foo, false);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'3foo'});
    tc.verifyEqual(out.d_3foo, true);
    tc.verifySize(args, [0 0]);
    
    % set
    opt.d_3foo = [];
    [out,args] = tb_optparse(opt, {});
    tc.verifyEqual(out.d_3foo, []);
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'3foo', [1 2 3]});
    tc.verifyEqual(out.d_3foo, [1 2 3]);
    tc.verifySize(args, [0 0]);
    
    % choose
    opt.d_3foo = {'this', 'that', 'other'};
    [out,args] = tb_optparse(opt, {});
    tc.verifyEqual(out.d_3foo, 'this');
    tc.verifySize(args, [0 0]);
    
    [out,args] = tb_optparse(opt, {'3foo', 'that'});
    tc.verifyEqual(out.d_3foo, 'that');
    tc.verifySize(args, [0 0]);
end