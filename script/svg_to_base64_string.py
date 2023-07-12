import base64
import os
import argparse
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outPath', type=str, help=
	'''input the outpath''',)
    args = parser.parse_args()
    return args.outPath
def png_to_base64(file=str()):
    #call get_args() to get outpath
    outpath = get_args() 
    outpath = outpath + "/report"
    filename = file.split('.')[0]
    format = file.split('.')[-1]
    file_path = outpath+"/"+file
    base64_path = outpath+'/base64'+'/'+filename+'.base64'
    if os.path.isfile(file_path):
        with open(file_path, "rb") as f:
            base64_data = base64.b64encode(f.read())
            s = base64_data.decode()
            base64_path_f = open(base64_path, 'w')
            if format =="svg":
                #base64_path_f.write('data:image/'+format+';base64,'+s)           
                base64_path_f.write(s) 
        
if __name__ == '__main__':
    path = get_args() 
    pictures = [] 
    for root,dir,files in os.walk(path):
        for f in files:
            if f.split('.')[-1] == "svg":
                pictures.append(f)
    #print(pictures)
    for p in pictures:
        png_to_base64(p) 
