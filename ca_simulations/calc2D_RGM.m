nSpec=6;

In=uint8(127*ones(Param.Array(1),Param.Array(2),1));

Init=zeros(Param.Array(1),Param.Array(2),nSpec);
Init(:,:,2)=Param.Array(3);
Init(:,:,3)=0;
Init(:,:,4)=Param.Array(4);
Init(:,:,5)=0;
Init(:,:,6)=Param.Array(7);
Init(:,:,1)=Param.Array(3)*Param.Array(6)*rand(Param.Array(1), Param.Array(2));
Init(:,:,2)=Param.Array(3)-Init(:,:,1);

[O,~]=maltesCA(Param,In,Init);

clear nSpec u