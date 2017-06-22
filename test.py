import numpy as np
from osgeo import ogr, gdal


def smoothingSwets(pixel, regWindow, avgWindow, nodata=None):
    #Check if entire pixel vector is nan, just return it
    if nodata and np.all(pixel==nodata):
        pixel.fill(np.nan)
        return pixel
    elif np.all(np.isnan(pixel)):
        return pixel
    
    #Get the shape of the original data
    originShape = pixel.shape
    
    #Reshape the data
    pixel = np.reshape(pixel, (len(pixel),1))
    pixel = pixel.astype(np.float32)
    
    #Replace no data values by nan if needed
    if nodata:
        pixel[pixel==nodata] = np.nan
    
    #Interpolate missing values if any
    #THIS SHOULD BE DONE BEFOREHAND FOR THE ENTIRE RASTER USING TEMPORAL AND SPATIAL INFERENCE
    #Get a mask of the nan values
    if np.isnan(pixel).any():
        pixel = np.reshape(pixel, len(pixel))
        nans = np.isnan(pixel)
        #Get the actual indices for these values and the non nan values
        notNans = ~nans
        notNans = notNans.nonzero()[0]
        nans = nans.nonzero()[0]
        #Interpolate
        pixel[nans]= np.interp(nans, notNans, pixel[notNans])
        #Reshape
        pixel = np.reshape(pixel, (len(pixel),1))
    
    #Create a vector of weights for the regression
    raw = pixel.ravel() #Flatten the array into a simple list
    raw = raw.tolist()
    weights = [0.5] #Consider the first point to be a local sloping point
    #Local peak = 1.5, local valley = 0.005, local sloping = 0.5
    for i in range(1, len(raw)-1, 1):
        if raw[i-1] < raw[i] and raw[i] > raw[i+1]:
            weights.append(1.5)
        elif ((raw[i-1] <= raw[i] and raw[i] <= raw[i+1]) or (raw[i-1] >= raw[i] and raw[i] >= raw[i+1])):
            weights.append(0.5)
        elif raw[i-1] > raw[i] and raw[i] < raw[i+1]:
            weights.append(0.005)
    weights.append(0.5) #Consider the last point to be a sloping point as well
    #Transform into an array
    weights = np.asarray(weights)
    weights = np.reshape(weights, pixel.shape)
    
    #Create a matrix with the data for this pixel and for the weights
    #For the data:
    #Each column will be the same
    dataRaw = np.repeat(pixel, len(pixel)+regWindow, axis=1)
    #Set to nan all data for each column besides the data used for the regressions
    ltri = np.tril_indices(n=len(pixel), k=-1, m=len(pixel)+regWindow) #Lower triangle indices below the diagonal
    dataRaw[ltri] = np.nan
    utri = np.triu_indices(n=len(pixel), k=0, m=len(pixel))
    dataRaw[:,(regWindow):][utri] = np.nan
    #Remove the first two and last 3 columns, since they don't have enough points for a regression
    dataRaw = dataRaw[:,2:(len(pixel)+regWindow-3)]
    
    #For the weights:
    weights = np.repeat(weights, len(pixel)+regWindow-5, axis=1)
    weights[np.isnan(dataRaw)] = np.nan
    
    #Create an empty array for the results of the regressions
    dataReg = np.zeros(dataRaw.shape)
    
    #Estimate the regression for each column of dataRaw
    for i in range(len(pixel)+regWindow-5):
        #Prepare regression data
        y = dataRaw[:,i][~np.isnan(dataRaw[:,i])] #dependent
        w = weights[:,i][~np.isnan(dataRaw[:,i])] #weights
        x = range(1,len(y)+1) #independent
        x = np.asarray(x)
        x = x.astype(np.float32)
        x = np.reshape(x, y.shape)
        #Estimate potential outliers
        out = np.divide(np.absolute(y-np.mean(y)),np.std(y)) #Will be an outlier if greater than 3 (only extreme outliers should be picked up)
        out = np.reshape(out, y.shape)
        #Remove outliers before regression
        yout = y[out<3]
        wout = w[out<3]
        xout = x[out<3]
        #Compute parameters of regression
        numerator = np.sum(wout)*np.sum(np.multiply(np.multiply(wout,xout),yout))-np.sum(np.multiply(wout,yout))*np.sum(np.multiply(wout,xout))
        denominator = np.sum(wout)*np.sum(np.multiply(wout,np.square(xout)))-np.square(np.sum(np.multiply(wout,xout)))
        b = np.divide(numerator, denominator)
        a = np.divide((np.sum(np.multiply(wout,yout))-b*np.sum(np.multiply(wout,xout))),np.sum(wout))
        #Compute the predicted values from the regression
        dataReg[:,i][~np.isnan(dataRaw[:,i])] = a + b*x
    
    dataReg[np.isnan(dataRaw)] = np.nan
    
    #Combination of the results from the regression for each point
    #Now we combine for each point the results from the moving regression window
    #We take into account the results from avg.window regression windows, centered around the point, unless we are at the edges
    smoothed = [] #Will hold the smoothed results
    t = int(np.floor(avgWindow/2.)) #number of predicted values from regressions to take into account around the center
    for i in range(len(pixel)):
        x = dataReg[i,:][~np.isnan(dataReg[i,:])]
        if i < np.floor(regWindow/2.) and len(x) < regWindow:
            center = int(np.floor(regWindow/2.))-(regWindow-len(x))
            res = np.mean(x[max(0,center-t):max(1,center+t+1)])
        else:
            res = np.mean(x[(int(np.floor(regWindow/2.))-t):(int(np.floor(regWindow/2.))+t+1)])    
        smoothed.append(res)
    
    smoothed = np.asarray(smoothed)
    smoothed = np.reshape(smoothed, originShape)
    
    return smoothed

#x = np.array([[9039,9021,8926,9009,8856,6887,8914,9076,8951,9015,8885,8877]])
x = np.array([[5740,5099,6533,6059,5595,6348,5771,3130,5685,6107,5893,6369]])
x = x.T
x = x.astype(np.float32)

print(smoothingSwets(x, 7, 3, nodata=-3000))

