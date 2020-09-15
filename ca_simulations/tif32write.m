function tif32write(imData,imName,imTags)
% A simple wrapper that saves 32-bit (single precision) floating values
% in a TIF-file. Usage:
%   tif32write(DATA, NAME, TAGS)
%     DATA -- N-dim array containing the image-Information. If N bigger
%             than 2, a stack of 2D-images will be saved
%     NAME -- char-string of the TIF-filename, old files will be 
%             overwritten
%   optional:
%     TAGS -- structure of valid TIFF-Tags, that are to be included in the 
%             TIF-file, but are not minimally necessary (e.g. XResolution).
%             Will be processed last, so use this to switch off or change
%             compression (default LZW at the moment)

    if (~isa(imData,'single')) imData=single(imData); end

    myTag.ImageLength = size(imData,1);
    myTag.ImageWidth = size(imData,2);
    if (ndims(imData)>3) imData=reshape(imData,myTag.ImageLength,myTag.ImageWidth,[]); end
    numPages =  size(imData,3);
    myTag.PageNumber=[1 numPages];

    myTag.Photometric = Tiff.Photometric.MinIsBlack;
    myTag.BitsPerSample = 32;
    myTag.SamplesPerPixel = 1;
    myTag.RowsPerStrip = 16;
    myTag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    myTag.Software = 'MATLAB';
    myTag.SampleFormat=Tiff.SampleFormat.IEEEFP;
    myTag.Compression=Tiff.Compression.LZW;

    if (nargin>=3 && isstruct(imTags))
        imTagNames=fieldnames(imTags);
        for u=1:numel(imTagNames)
            myTag=setfield(myTag,imTagNames{u},getfield(imTags,imTagNames{u}));
        end
    end
    
    t=Tiff(imName,'w');
    t.setTag(myTag);
    t.write(imData(:,:,1));
    
    for u=2:numPages;
        t.writeDirectory();
        myTag.PageNumber=[u numPages];
        t.setTag(myTag);
        t.write(imData(:,:,u));
    end
    
    t.close();
end