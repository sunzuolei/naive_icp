function [TimeLaser, Range, Intensity] =ReadLaserData(file)

    load(file);
    TimeLaser=double(TLsr)/1000;
    
    Mask13 = uint16(213 -1) ;                              %Masks to get the range and intensity vectors
    MaskA  = bitcmp(Mask13,16) ;
    Range=zeros(size(LASER)); Intensity=zeros(size(LASER));
    for i=1:size(LASER,1)
        Laser1=LASER(i,:);
        RR = double(  bitand( Mask13,Laser1) ) ;             %range
        II  = uint16(  bitand( MaskA ,Laser1) ) ;            %intensity, >0, high intensity
        RR = RR/100 ;           %cm ---> metros
        Range(i,:)=RR;
        Intensity(i,:)=II;
    end
return