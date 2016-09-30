function X=Contract(tensors,contractions)

    numcont=max(cell2mat(contractions));

    %fprintf('numcont is %d\n', numcont);

    table=zeros(numcont,2);
    for k=1:length(contractions)
        c=contractions{k};
        for i=c(c>0)
            if table(i,1)==0
                table(i,1)=k;
            else
                table(i,2)=k;
            end
        end
    end
    
    contlist=1:numcont;

    while ~isempty(contlist);
        indX=table(contlist(1),1);
        cX=contractions{indX};
        X=tensors{indX};
        tensors{indX}=[];
        dimX=size(X);
        if length(cX)>length(dimX)
            dimX=[dimX ones(1,length(cX)-length(dimX))];
        end
        indY=table(contlist(1),2);
        cY=contractions{indY};
        Y=tensors{indY};
        tensors{indY}=[];
        contractions{indY}=[];
        dimY=size(Y);
        if length(cY)>length(dimY)
            dimY=[dimY ones(1,length(cY)-length(dimY))];
        end

        [contrlines poscX poscY]=findcontr(cX,cY);
        for ic=contrlines
            contlist(contlist==ic)=[];
        end
        posrX=1:length(cX);posrX(poscX)=[];
        posrY=1:length(cY);posrY(poscY)=[];
        if length(cX)>1
            X=permute(X,[posrX poscX]);
        end
        X=reshape(X,[prod(dimX(posrX)) prod(dimX(poscX))]);
        if length(cY)>1
            Y=permute(Y,[poscY posrY]);
        end
        Y=reshape(Y,[prod(dimY(poscY)) prod(dimY(posrY))]);
        X=X*Y;
        clear Y;
        dimX=[dimX(posrX) dimY(posrY)];
        cX=[cX(posrX) cY(posrY)];
        if isempty(dimX)
            dimX=[1 1];
        elseif length(dimX)==1
            dimX=[dimX 1];
        end
        tensors{indX}=reshape(X,dimX);
        contractions{indX}=cX;
        table(table==indY)=indX;
    end

    X=[];
    for indX=1:length(tensors)
        if ~isempty(tensors{indX})
            if isempty(X)
                X=tensors{indX};
                cX=contractions{indX};
                dimX=size(X);
                if length(cX)>length(dimX)
                    dimX=[dimX ones(1,length(cX)-length(dimX))];
                end
            else
                Y=tensors{indX};
                cY=contractions{indX};
                dimY=size(Y);
                if length(cY)>length(dimY)
                    dimY=[dimY ones(1,length(cY)-length(dimY))];
                end
                X=reshape(X,[numel(X) 1]);
                Y=reshape(Y,[1 numel(Y)]);
                dimX=[dimX dimY];
                cX=[cX cY];
                X=reshape(X*Y,dimX);
                clear Y;
            end
        end
    end

    if length(cX)>1
        [indX,orderX]=sort(-cX);
    else
        orderX=[1 2];
    end
    X=permute(X,orderX);

end


function [contrlines poscX poscY]=findcontr(cX,cY)

    X1=cX.'*ones(1,length(cY));
    X2=ones(length(cX),1)*cY;

    ind=find(X1==X2)';
    contrlines=X1(ind);
    poscX=mod(ind-1,length(cX))+1;
    poscY=floor((ind-1)/length(cX))+1;

end