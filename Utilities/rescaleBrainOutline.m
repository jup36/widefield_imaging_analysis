function rescaledMap = rescaleBrainOutline(origMap, rescaledSize)
    rescaledMap = imresize(origMap, rescaledSize, 'nearest'); 
end


