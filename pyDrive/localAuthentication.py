# import dependencies

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
import sys

# instantiate GoogleAuthentication Object

GAuth = GoogleAuth()

# prompt user for name

name = input("Please type your name: ")

# Use your local machines web server at localhost port 8080 to authenticate

GAuth.LocalWebserverAuth()

# save authentication code in a credentials file to be used on the server

GAuth.SaveCredentialsFile(name+"sCredentials.txt")

