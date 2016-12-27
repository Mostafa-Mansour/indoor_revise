function [ wallNum ] = getWallNum( x1,x2,k,D )

if k>360
    k=k-360;
end

switch k
    case 0
        wallNum=3;
    case 90
        wallNum=4;
    case 180
        wallNum=1;
    case 270
       wallNum=1;
    case 360
        wallNum=3;
    otherwise
        
        if k>0&&k<90
            m=tand(90-k);
            if m<(D(2)-x2)/(D(1)-x1)
                wallNum=4;
            else
                wallNum=3;
            end
        
    
        elseif k>90 && k<180
            m=tand(180-(k-90));
            if m<-x2/(D(1)-x1)
                wallNum=1;
            else
                wallNum=4;
            end
                
        
        elseif k>180 && k<270
            m=tand(90+(360-k));
            if m>x2/x1
                wallNum=1;
            else
                wallNum=2;
            end

                
        else
            m=tand(90+(360-k));
            if m>(D(2)-x2)/(-x1)
                wallNum=2;
            else
                wallNum=3;
            end          
               
       
        
       end
    
   
end  


end

