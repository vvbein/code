function [error] = node_optimization(x,P,T,inputps,outputps)
    numpen0 = round(x_opt(1));
    numpen1 = round(x_opt(2));

    %% unsupervised pre-training
    fprintf(1,'Pretraining Layer 1 with RBM: %d-%d ',numdims,numpen0);
    numhid=numpen0;
    restart=1;
    rbm1;
    vishid1=vishid;hidrecbiases=hidbiases;

    fprintf(1,'\nPretraining Layer 2 with RBM: %d-%d ',numhid,numpen1);
    batchdata=batchposhidprobs;
    numhid=numpen1;
    restart=1;
    rbm1;
    hidpen=vishid; penrecbiases=hidbiases; hidgenbiases=visbiases;
    %trained RBMs are used for the initialized the DBN weights
    w1=[vishid1; hidrecbiases];
    w2=[hidpen; penrecbiases];
    %% supervised regression
    N1 = size(P,1);
    digitdata = [P ones(N1,1)];

    w1probs = 1./(1 + exp(-digitdata*w1));
    w1probs = [w1probs  ones(N1,1)];

    w2probs = 1./(1 + exp(-w1probs*w2));
    w2probs = [w2probs  ones(N1,1)];

    H= w2probs'; 
    %% BP fine-tunning
    min_v= min(numpen0, numpen1);
    inputn=H(1:min_v,1:10:1800);
    outputn=T';
    net=newff(inputn,outputn,[numpen0,numpen1]);
    net.trainParam.epochs=100;
    net.trainParam.lr=0.1;
    net.trainParam.goal=0.03;
    net.trainParam.max_fail = 200;
    net.trainParam.showWindow = false;
    net.trainParam.showCommandLine = false;
    net=train(net,inputn,outputn);
    input_test1=data_none(1801:10:2700,1:8);     
    output_test1=data_none(1801:10:2700,9);
    input_test=input_test1';
    output_test=output_test1';
  
    test_x=mapminmax('apply',input_test,inputps,0,1)';

    N1 = size(test_x,1);
    digitdata1 = [test_x ones(N1,1)];

    w1probs1 = 1./(1 + exp(-digitdata1*w1));
    w1probs1 = [w1probs1  ones(N1,1)];

    w2probs1 = 1./(1 + exp(-w1probs1*w2));

    test_x1 = w2probs1';
    an=sim(net,test_x1);

BPoutput=mapminmax('reverse',an,outputps);
error=BPoutput-output_test;
error=sum(error)./90*100
end