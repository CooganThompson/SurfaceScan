function plotbrokenlines(XCoords,YCoords,varargin)

%for i=1:length(XCoords)
    
    Ydis=abs(diff(YCoords));
    skipped=[0 (Ydis>100)'];
    blocks = cumsum(skipped);
    if sum(blocks)>0
        
        indexofchange=diff(blocks);
        point=find(indexofchange==1);
        indexofchange(point+1:end+1)=indexofchange(point:end);
        
        IndexofFirst=find(indexofchange==1,1);
        FirstBeta=YCoords(IndexofFirst);
        SecondBeta=YCoords(IndexofFirst+1);
        
        if FirstBeta<180
            FirstBeta=FirstBeta+360;
            flag=1;
        elseif SecondBeta<180
            SecondBeta=SecondBeta+360;
            flag=2;
        else
            error('Something is wrong with the points')
        end
        
        Yvalues=XCoords(indexofchange==1);
        slope=diff(Yvalues)/(SecondBeta-FirstBeta);
        
        if flag==1
            NewYvalue=slope*(360-FirstBeta)+Yvalues(1);
        elseif flag==2
            NewYvalue=slope*(360-SecondBeta)+Yvalues(2);
        end
        
        XCoords(point+3:end+2)=XCoords(point+1:end);
        YCoords(point+3:end+2)=YCoords(point+1:end);
        
        XCoords(point+1:point+2)=NewYvalue;
        
        if flag==1
            YCoords(point+1:point+2)=[0 360];
        elseif flag==2
            YCoords(point+1:point+2)=[360 0];
        end
        
        Ydis=abs(diff(YCoords));
        skipped=[0 (Ydis>100)'];
        blocks = cumsum(skipped);
        
    end
    
    for j=0:blocks(end)
        % plot(XCoords{i}(blocks==j),YCoords{i}(blocks==j),...
        %   'linestyle','-','color',Color,'linewidth',2);
        plot(XCoords(blocks==j),YCoords(blocks==j),varargin{:});
        %                    hold on
    end
end