#from sys import *
import scipy as sc
import numpy as np
#from matplotlib import *
import matplotlib.pyplot as plt
from scipy import io
from scipy.io import wavfile
from scipy import ndimage
from pylab import specgram, psd
from matplotlib import mlab
from itertools import *
from scipy import signal
from array import array


def impwav(a):
 """Imports a wave file as an array where a[1] 
 is the sampling frequency and a[0] is the data"""
 out=[]
 wav = sc.io.wavfile.read(a)
 out=[wav[1],wav[0]]
 return out

'''def impcbin(a):
 """Imports a cbin as an array where a[1] is the sampling freq and a[0] is the data"""
 file=open(a,'r')
 file=file.read()
 a=array('h')
 a.fromstring(file)
 return[np.asarray(a),32000]'''

def impcbin(a):
 """Imports a cbin as an array where a[1] is the sampling freq and a[0] is the data"""
 file=open(a,'r')
 file=file.read()
 a=np.fromstring(file[1:-1],dtype=np.int16)
 out=[a,32000]
 return out

def outwave(filename,array):
 """Exports a numpy array (or just a regular python array) 
 as a wave file. It expects an array of the following format: (sample np.sqrt(sum((np.array(a)-np.array(b))**2))freq, data)"""
 sc.io.wavfile.write(filename,array[0],array[1])

def spec1d(a):
    "plots a spectrogram of a"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.specgram(a,Fs=32000)
    #if a[1]==None: ax.specgram(a[0])
    #else: ax.specgram(a[0],Fs=a[1])
    plt.show()

#def ryth1d(a):
#    "plots a spectrogram of a.  4_13_11 I actually have no idea what I was tyring to do here."
#    fig=plt.figure()
#    ax=fig.add_subplot(111)np.sqrt(sum((np.array(a)-np.array(b))**2))
#    ax.specgram(a,NFFT=(len(a)/2),Fs=32000)
#    #if a[1]==None: ax.specgram(a[0])
#    #else: ax.specgram(a[0],Fs=a[1])
#    plt.show()

def plot1d(a):
    "plots a"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(a)
    plt.show()

def scat1d(a):
    "scatterplots a. Expects a to be a ziped string of the two values"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.scatter(a)
    plt.show()

def hist(a,bins):
    "plots a histogram of a"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.hist(a,bins)
    plt.show()

def histlog(a,bins):
    "plots a histogram of a"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.hist(a,bins,log=True)
    plt.show()

def threshold(a,thresh=None):
 """Returns a thresholded array of the same length as input
 with everything below a specific threshold set to 0.
 By default threshold is sigma."""
 if thresh==None: thresh = sc.std(a)
 out=np.where(abs(a)>thresh,a,np.zeros(a.shape))
 return out

def rthreshold(a, thresh=None):
 """Returns a reverse thresholded array of the same length as input
 with everything above a specific threshold set to 0.
 By default threshold is sigma."""
 if thresh==None: thresh = sc.std(a)
 out=np.where(abs(a)>thresh,np.zeros(a.shape),a)
 return out

def mask(a,thresh=None):
 """Returns a masnp.sqrt(sum((np.array(a)-np.array(b))**2))k array 
 of the same length as input
 with everything below a specific threshold set to 0 and
 everything above that threshold set to 1.
 By default threshold is sigma."""
 freq = a[1]
 a = a[0]
 if thresh==None: thresh = 5*np.std(a)
 out=np.where(abs(a)>thresh,np.ones(a.shape),np.zeros(a.shape))
 return out


def weinent(a,range=None,nfft=None,fs=None):
 "returns the weiner entropy of an array"
 if range==None: range=(500,8000)
 if nfft==None: nfft=256
 if fs==None: fs=32000
 segstart=int(range[0]/(fs/nfft))
 segend=int(range[1]/(fs/nfft))
 a=np.array(a)
 a=mlab.psd(a,Fs=fs,NFFT=nfft)
 a=np.array(a[0][segstart:segend])
 log_a=np.log(a)
 out=((np.exp(log_a.mean()))/(a.mean()))
 return out

#def window(seq, n):
#    "Returns a string of windows (of width n). 4_13_11 getting rid of this because it is a sloppy solution to this problem."
#    it = iter(seq)
#    out = tuple(islice(it, n))
#    if len(out) == n:
#     return out

def mfreqz(b,a=1):
 """Plots the frequency and phase response of a filter.""" 
 w,h = sc.signal.freqz(b,a)
 fig=plt.figure()
 h_dB = 20 * np.log10 (abs(h))
 ax=fig.add_subplot(211)
 ax.plot(w/max(w),h_dB,'k')
 ax.set_ylim(-150, 5)
 ax.set_ylabel('Magnitude (db)')
 ax.set_xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
 ax.set_title(r'Frequency response')
 ay=fig.add_subplot(212)
 h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h)))
 ay.plot(w/max(w),h_Phase,'k')
 ay.set_ylabel('Phase (radians)')
 ay.set_xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
 ay.set_title(r'Phase response')
 fig.subplots_adjust(hspace=0.5)
 plt.show()

def filtercall(a):
 b=sc.signal.iirdesign([0.5,0.7], [0.46,0.75], 1, 60, ftype='ellip')
 out=sc.signal.filtfilt(b[0],b[1],a)
 return (out)

def filtersongpb(a):
 out=[]
 b=sc.signal.iirdesign(wp=[0.04,0.63], ws=[0.02,0.67], gpass=1, gstop=60,ftype='ellip')
 out.append(sc.signal.filtfilt(b[0],b[1],a[0]))
 out.append(a[1])
 return (out)

def filtersong(a):
 out=[]
 b=sc.signal.iirdesign(wp=0.04, ws=0.02, gpass=1, gstop=60,ftype='ellip')
 out.append(sc.signal.filtfilt(b[0],b[1],a[0]))
 out.append(a[1])
 return (out)

def smoothrect(a,window=None,freq=None):
 """smooths and rectifies a song.  Expects a file format (data,samplerate).
    If you don't enter a smoothing window size it will use 2ms as a default."""
 if freq== None:
  #freq=32000
  freq=a[1]
 if window == None: window=2
 le=round(freq*window/1000)
 h=np.ones(le)/le
 smooth= np.convolve(h,abs(a[0]))
 offset = round((len(smooth)-len(a[0]))/2)
 smooth=smooth[(1+offset):(len(a[0])+offset)]
 return (smooth, freq)

def comb (a,size):
 """FIX ME! too slow!this is intended to give you a mask that gives 
 you larger objects it expects a binary mask as input.  It is 
 equivalent to a convolution of a mask and then thresholding it. 
 the expected size is in samples"""
 return [1 if max(a[n:n+size])==1 else 0 for n in range(len(a)-size)]

def remobjmsk(msk,size=None,freq=None):
 """FIX ME! too slow! removes objects that are too short from a mask"""
 if size== None: size=30
 if freq== None: freq=32000
 sizeinsamps=round(freq*size/1000)
 label=(sc.ndimage.label(msk))
 label=[list(label[0]),label[1]]
 keepsyls=[x for x in range(label[1]) if label[0].count(x)>sizeinsamps]
 outmsk=[x if x in keepsyls else 0 for x in label[0]]
 return outmsk 

def smooth(a,window=None,freq=None):
 """smooths a song.  Expects a file format (data,samplerate).
    If you don't enter a smoothing window size it will use 2ms as a default."""
 if window == None: window=2
 if freq== None: freq=32000
 le=round(freq*window/1000)
 h=np.ones(le)/le
 smooth= np.convolve(h,a)
 offset = round((len(smooth)-len(a))/2)
 smooth=smooth[(1+offset):(len(a)+offset)]
 return smooth

def specplot(a):
 fig=plt.figure()
 ay=fig.add_subplot(111)
 #ay.specgram(a[0],Fs=a[1])
 ay.plot(smoothrect(filtersong(a),20), 'k', (mask(smoothrect(filtersong(a),20),2*(np.abs(np.median(a[0]))))),'b')
 plt.show()

def songseg(a, sigma=None):
 """returns a mask of an array where everything below 2sigma is zero and anything above 2
 sigma is a number 1 through N which is the number of the syllable.
 This is really helpful for segmenting data. This expects the song to have been imported using
 impwav."""
 if sigma==None: sigma=2
 label=sc.ndimage.label(mask(smoothrect(filtersong(a),20),sigma*(np.abs(np.median(a[0])))))
 return label

def songsegrange(a, cutoff=None,smwindow=None):
 """returns a mask of an array where everything below 2sigma is zero and anything above 2
 sigma is a number 1 through N which is the number of the syllable.
 This is really helpful for segmenting data. This expects the song to have been imported using
 impwav."""
 if smwindow==None: smwindow=20
 if cutoff==None: cutoff=2
 label=sc.ndimage.label(mask(smoothrect(filtersong(a)[0],smwindow),cutoff))
 return label

def sylablelen(a):
 """returns an array with the length of all sylables in wave file a.  Expects a
 to have been created using impwav.  At some point I should probably rewrite this to be faster."""
 lab=songseg(a)
 freq=a[1]
 sylno=lab[1]
 inc=1
 out=[]
 lst=list(lab[0])
 while inc<=sylno:
  len=lst.count(inc)
  out.append(len)
  inc=inc+1
 out=out/freq
 return out

def intersyllen(a):
 """returns an array with all intersylable lengths in wave file a.  Expects a
 to have been created using impwav."""
 msk=mask(smoothrect(filtersong(a),20),2*np.median(a[0]))
 freq=a[1]
 lst=list(msk)
 lst.reverse()
 endind=lst.index(1)
 lst.reverse()
 lst=lst[lst.index(1):(len(lst)-(endind))]
 lst=np.array(lst)
 lst=np.where(lst==0,np.ones(len(lst)),np.zeros(len(lst)))
 lst=sc.ndimage.label(lst)
 out=[]
 sylno=lst[1]
 inc=1
 lst=list(lst[0])
 while inc<=sylno:
  leng=lst.count(inc)
  out.append(leng)
  inc=inc+1
 out=np.float32(out)
 out=out/(int(freq)/1000)
 return out

def sylrate(a):
 """Returns the rate of sylables in a given song in sylables/sec"""
 lab=songseg(a)
 msk=mask(smoothrect(filtersong(a),20),2*np.median(a[0]))
 lst=list(msk)
 lst.reverse()
 endind=lst.index(1)
 lst.reverse()
 lst=lst[lst.index(1):(len(lst)-(endind))]
 songleng= np.float32(len(lst))/np.float32(a[1])
 out=lab[1]/songleng
 return out

def songleng(a):
 """Returns the song length in sylables"""
 lab=songseg(a)
 return lab[1]

def ffcalc(a, freq=None):
 """Returns the fundamental frequency of a string, a. Expects raw data,  the default frequency is 32,000. This uses brute force correlation which is slow for large data sets but more accurate than fft based methods. Returns the data in wavelength"""
 if freq==None: freq=32000
 corr=sc.correlate(a,a,mode='same')
 corr=corr[(len(corr)/2):(len(corr)-len(corr)/4)]
 dat=np.diff(np.where(np.diff(corr)>0,1,0))
 out=float(freq)/float(((list(dat)).index(-1)))
 return out

def ffcalcfft(a, freq=None):
 """Returns the fundamental frequency of a string, a. Expects raw data,  the default frequency is 
 32000. This method uses ffts."""
 if freq==None: freq=32000
 fft=sc.fftpack.fft(a)
 corr=sc.fftpack.ifft(fft*fft.conjugate())
 corr=corr[:(len(corr)/4)]
 dfff=np.diff(corr)
 dat=np.diff(np.where(dfff>0,1,0))
 if -1 not in list(dat): out=freq
 else:
  first=(((list(dat)).index(-1)))
  slope=(dfff[(first+1)]-dfff[first])
  out=(slope*first-dfff[first])/slope
 out=freq/out
 return out

def ffqualitymask(a, freq=None):
 """NOT DONE! Returns the fundamental frequency of a string, a. Expects raw data,  the default frequency is 
 32000. This method uses ffts."""
 if freq==None: freq=32000
 fft=sc.fftpack.fft(a)
 corr=sc.fftpack.ifft(fft*fft.conjugate())
 corr=corr[:(len(corr)/4)]
 dfff=np.diff(corr)
 dat=np.diff(np.where(dfff>0,1,0))
 if -1 not in list(dat): out=freq
 else:
  first=(((list(dat)).index(-1)))
  slope=(dfff[(first+1)]-dfff[first])
  out=(slope*first-dfff[first])/slope
 out=freq/out
 return out

def ffcalcfftqual(a, freq=None):
 """Returns the fundamental frequency and the amplitude of the autocorrelation at the first peak, of a string, a.
 Expects raw data,  the default frequency is 32000. This method uses ffts."""
 if freq==None: freq=32000
 fft=sc.fftpack.fft(a)
 corr=sc.fftpack.ifft(fft*fft.conjugate())
 corr=corr[:(len(corr)/4)]
 dfff=np.diff(corr)
 dat=np.diff(np.where(dfff>0,1,0))
 if -1 not in list(dat): out=(float(0),float(0))
 elif list(dat).index(-1)<2: out=(float(0),float(0))
 else:
  first=(((list(dat)).index(-1)))
  slope=(dfff[(first+1)]-dfff[first])
  out=(slope*first-dfff[first])/slope
  out=(freq/out,corr[first])
 return out

def ffcalcfftqual2(a, freq=None):
 """Returns the fundamental frequency of a string, and the amplitude of the correlation at the best peak (of the first 4).
 the default frequency is 32000. This method uses ffts."""
 if freq==None: freq=32000
 fft=sc.fftpack.fft(a)
 corr=sc.fftpack.ifft(fft*fft.conjugate())
 corr=corr[:(len(corr)/4)]
 dfff=np.diff(corr)
 dat=np.diff(np.where(dfff>0,1,0))
 if -1 not in list(dat): out=(float(freq),float(0))
 else:
  dat=np.where(dat<0,1,0)
  dat=list(dat) 
  (lab, inc)=sc.ndimage.label(dat)
  if inc >4: inc=4
  arr=[]
  while inc!=0:
   first=[]
   pos=((list(lab)).index(inc))-1
   first.append(pos)
   first.append(corr[pos])
   arr.append(first)
   inc=inc-1
  arr.sort(key=lambda x:x[1])
  first=arr[len(arr)-1][0]
  slope=(dfff[(first+1)]-dfff[first])
  out=(slope*first-dfff[first])/slope
  out=freq/out
  out=(out,(arr[len(arr)-1][1]))
 return out

def ffcalcfft2(a, freq=None):
 """Returns the fundamental frequency of a string, corrected to get the best correlation not the first peak
 Default freq=32000. This method uses ffts."""
 if freq==None: freq=32000
 fft=sc.fftpack.fft(a)
 corr=sc.fftpack.ifft(fft*fft.conjugate())
 corr=corr[:(len(corr)/4)]
 dfff=np.diff(corr)
 dat=np.diff(np.where(dfff>0,1,0))
 if -1 not in list(dat): out=freq
 else:
  dat=np.where(dat<0,1,0)
  dat=list(dat)
  (lab, inc)=sc.ndimage.label(dat)
  if inc >4: inc=4
  arr=[]
  while inc!=0:
   first=[]
   pos=((list(lab)).index(inc))-1
   first.append(pos)
   first.append(corr[pos])
   arr.append(first)
   inc=inc-1
  arr.sort(key=lambda x:x[1])
  first=arr[len(arr)-1][0]
  slope=(dfff[(first+1)]-dfff[first])
  out=(slope*first-dfff[first])/slope
  out=freq/out
 return out

def ffprofilequal(a, freq=None, window=None):
 """returns a string of local estimates of the fundamental frequency of a string.
    It also returns the amplitude of the highest peak in the autocorrelation as a quality measure.
    By default, freq=32000 and window=256."""
 if freq==None: freq=32000
 if window==None: window=256 
 out = [ffcalcfftqual(a[x:x+(window-1)],freq) for x in range(len(a))] 
 return out

def ffprofilequalmask(a, freq=None, window=None):
 """returns a string of local estimates of the fundamental frequency of a string.
    It also returns the amplitude of the highest peak in the autocorrelation as a quality measure.
    By default, freq=32000 and window=256."""
 if freq==None: freq=32000
 if window==None: window=256
 out = [ffcalcfftqual(a[x:x+(window-1)],freq) for x in range(len(a))]
 dat,sigma=zip(*out)
 sigma=sc.std(sigma)
 out=[out[0] for x in out if out[1] >sigma]
 return out

def ffprofilequal2(a, freq=None, window=None):
 """returns a string of local estimates of the fundamental frequency of a string correcting for harmonic errors.
    It also returns the amplitude of the highest peak in the autocorrelation as a quality measure.
    By default, freq=32000 and window=256."""
 if freq==None: freq=32000
 if window==None: window=256
 out = [ffcalcfftqual2(a[x:x+(window-1)],freq) for x in range(len(a))]
 return out

def ffprofile2(a, freq=None, window=None):
 """returns a string of local estimates of the fundamental frequency of a string correcting for harmonics
    By default, freq=32000 and window=256."""
 if freq==None: freq=32000
 if window==None: window=256
 out = [ffcalcfft2(a[x:x+(window-1)],freq) for x in range(len(a))]
 return out

def ffprofilefft(a, freq=None, window=None):
 """Returns a string of local estimates of the fundamental frequency of a string."""
 if freq==None: freq=32000
 if window==None: window=256
 out = [ffcalcfft(a[x:x+(window-1)],freq) for x in range(len(a))]
 return out

#def ffxcorr(a,b):
# """TO DO. Returns the cross correlation of a and b.  Expects raw data"""

def ffacorr(a):
 """Returns the autocorrelation of a. Expects raw data"""
 z=np.zeros(2*len(a))
 z[:len(a)]=a
 fft=sc.fft(z)
 out=sc.ifft(fft*sc.conj(fft))
 return (out[:len(out)/2])

def norm(a):
 """normalizes a string by it's average"""
 a=(np.array(a)-np.average(a))/np.std(a)
 return a

def maxnorm(a):
 """normalizes a string by it's max"""
# a=np.array(a)
 a=a/a.max()
 return a

def norment(a):
 """returns an entropy calculation for an array normalized to between 0 and 1"""
 a=a/np.average(a)
 out=-100*sum(a*np.log(a))/np.log(len(a))
 return out

def normmaxent(a):
 """returns an entropy calculation for an array normalized to between 0 and 1"""
 a=a/max(a)
 out=-100*sum(a*np.log(a))/np.log(len(a))
 return out

def entropy(a):
 """retunrs an entropy calculation for an array"""
 out=-100*sum(a*np.log(a))/np.log(len(a))
 return out

def shanonent(x):
 '''returns shannon entropy'''
 s=set(x)
 ents=[]
 for y in s:
  prb=x.count(y)/float(len(x))
  ents.append(prb*np.log(prb))
 return(-sum(ents))

def vshannonent(x):
 '''vectorized shannon entropy'''
 s=set(x)
 b=float(len(s))
 prbs=np.array([x.count(y)/float(len(x)) for y in s])
 H=-sum(prbs*(np.log(prbs)/np.log(b)))
 return H

def prb(y,x):
 '''returns probability of y in vector x'''
 return (list(x).count(y)/float(len(x)))

def jointent(y,x):
 '''returns conditional shannon entropy of y given x'''
 s=set(x+y)
 b=float(len(s))
 xprbs=np.array([x.count(n)/float(len(x)) for n in s])
 jprbs=np.array([prb(m,[y[nn] for nn in range(len(x)) if x[nn]==n])*prb(n,x) for n in s for m in s])
 jprbs=[x for x in jprbs if x!=0]
 HXY=-sum(jprbs*(np.log(jprbs)/np.log(b)))
 HX=-sum(xprbs*(np.log(xprbs)/np.log(b)))
 return HXY-HX


def window(func,a, window=None):
 """produces a string which is the application of func on sliding window on string a."""
 if window==None: window=256
 out = [func(a[x:x+(window-1)]) for x in range(len(a))]
 return out

def songfinder(song):
 """identifies an imported wav file as a song.  Returns True if it is a song
 and False if it isn't.  Earliest version based on mean amplitude"""
 powerspecd=plt.psd(filtersong(song)[0])
 ent=np.mean((powerspecd[0][10:100]))
 return ent>100000

def songfinder2(song):
 psd = signal.welch(filtersong(song)[0], fs=song[1])
 psd_var = np.log10(np.var(psd[1][27:65])) # equivalent to ~5K - 10K
 ent = weinent(song[0], fs=song[1])
 hifreq = np.log10(np.mean(psd[1][110:]))
 mean_pow = np.log10(np.mean(psd))
 return psd_var > 27 and hifreq > 9 and ent > .5 and mean_pow > 14

def eucliddist(a,b):
 """calculates a euclidean distance between two equal length arrays"""
 return np.sqrt(sum((np.array(a)-np.array(b))**2))

def sqformdistmat(array):
 """creates a squareform distance matrix for an array.  Currently only
 uses euclidean dist"""
 out=[]
 for x in array:
  intout=[]
  for y in array:
    intout.append(eucliddist(x,y))
  out.append(intout)
 return out

def findobject(file):
 value=(np.average(file))
 thresh=threshold(file,value)
 label=(sc.ndimage.label(thresh)[0])
 objs=sc.ndimage.find_objects(label)
 return(objs,thresh)

def psd_fft(a,fs=None):
 if fs==None:fs=32000
 n=float(len(a))
 ffta=sc.fft(a)
 ffta=ffta[1:n/2+1]
 ffta=np.array(ffta)
 psda=(1/(fs*n))*abs(ffta)
 psda[2:-1]=2*psda[2:-1]
 return(psda)

def psdrange_fft(a,fs=None,frange=None):
 if frange==None:frange=(650,10000)
 if fs==None:fs=32000
 n=float(len(a))
 ffta=sc.fft(a)
 ffta=ffta[1:n/2+1]
 ffta=np.array(ffta)
 psda=(1/(fs*n))*abs(ffta)
 psda[2:-1]=2*psda[2:-1]
 freq=np.arange(frange[0],frange[1],fs/n)
 psda=psda[frange[0]/(fs/n)-1:frange[1]/(fs/n)]
 return(psda)

def psd_lomb_scargle(a,fs=None,frange=None):
 if frange==None:frange=(650,10000)
 if fs==None:fs=32000
 n=float(len(a))
 x=np.linspace(0,1/32,n)
 f=np.linspace(500,10000,10)
 psd=sc.signal.lombscargle(x,a[0],f)


def pltpsd(a,fs=None,frange=None):
 if frange==None:frange=(650,10000)
 if fs==None:fs=32000
 n=float(len(a))
 ffta=sc.fft(a)
 ffta=ffta[1:n/2+1]
 ffta=np.array(ffta)
 psda=(1/(fs*n))*abs(ffta)
 psda[2:-1]=2*psda[2:-1]
 #freq=np.arange(0,(fs/2),fs/n)
 freq=np.arange(frange[0],frange[1],fs/n)
 psda=psda[frange[0]/(fs/n)-2:frange[1]/(fs/n)]
 npsda=norm(psda)
 smnpsda=smooth(npsda,20)
 print len(smnpsda)
 print len(freq)
 print len(psda)
 plt.plot(freq,(smnpsda))
 plt.title('psd using fft')
 plt.xlabel('freq (Hz)')
 plt.ylabel('power(dB)')
 plt.show()
 return(psda)

def it_otsu(input):
 input=np.array(input)
 currsong=input
 thrs=[]
 for x in range(10):
  currwsig=100
  currthr=100
  for thr in np.array(range(100))/1000.0:
   above=[x for x in (currsong>thr)*currsong if x>0]
   below=[x for x in (currsong<thr)*currsong if x>0]
   alen=float(len(above))
   blen=float(len(below))
   tlen=float(len(input))
   wsig=(alen/tlen*np.var(above))+(blen/tlen*np.var(below))
   #wsigs.append(wsig)
   #thrs.append(thr)
   if wsig<currwsig:
    currwsig=wsig
    currthr=thr
  thrs.append(currthr)
  currsong=[x if x < currthr else currthr for x in currsong]
 return min(thrs)


def songfinder(song):
 """identifies an imported wav file as a song.  Returns True if it is a song
 and False if it isn't.  Earliest version based on mean amplitude"""
 powerspecd=plt.psd(filtersong(song)[0])
 ent=np.mean(powerspecd[0][10:100])
 return ent>100000

def sylrate(a,thr):
 """Returns the rate of sylables in a given song in sylables/sec"""
 fsong,frq=filtersong(a)
 smfilt=((smoothrect(smoothrect(fsong,2,freq=frq),12,freq=frq)))
 #smfilt=((smoothrect(fsong,10,freq=frq)))
 #thr=2*max(smfilt[500:1500])
 #thr=2*it_otsu(smfilt[6400:128000:1000])
 msk=mask(smfilt,thr)
 #msk=mask(smfilt,0.275*np.std(smfilt))
 lab=sc.ndimage.label(msk)
 objs=sc.ndimage.find_objects(lab[0])
 objs=[x for x in objs if len(a[0][x[0]])>(5*frq/1000.0)]
 onsets=[x[0].start for x in objs]
 offsets=[x[0].stop for x in objs]
 if len(offsets)>0 and len(onsets)>0:
  songleng=(offsets[-1]-onsets[0])/float(frq)
 else:
  return 'nodata'
 syllens=np.array(offsets)-np.array(onsets)
 isis=np.array(onsets)[1:]-np.array(offsets)[:-1]
 syllens=syllens/float(frq) #added may 2015
 isis=isis/float(frq) #added may 2015
 sylpersec=len(objs)/songleng
 outobjs=zip(onsets,offsets)
 return sylpersec,syllens,isis,outobjs
