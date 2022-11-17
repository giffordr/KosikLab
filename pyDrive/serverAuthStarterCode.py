# import dependencies

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive


# instantiate GoogleAuthentication Object

GAuth = GoogleAuth()

# prompt user for name

name = input("Please type your name: ")

# load Credentials and Authorize Google Drive Connection

GAuth.LoadCredentialsFile(name+"sCredentials.txt")
Gauth.Authorize()

# instantiate a GoogleDrive object

Drive = GoogleDrive(GAuth)

# test sending file to your google drive

File1 = Drive.CreateFile({'title': 'Hello.txt'})
File1.SetContentString('Hello World!')
File1.Upload()



