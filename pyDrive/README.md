# KosikLab

This folder can be used to help mount a Google Drive to the CNSI server, to work with files from, or export to, your Google Drive.

We will use pydrive and Google Drive API to make this happen. You will first want to make sure that you have pydrive downloaded. Pull the Kosik Lab git to you local machine. Change directory to this folder in the Kosik Lab repo "pyDrive". Install pydrive using pip:

``` python -m pip install pydrive ```

Check to make sure that pydrive is installed correctly. One way to do that is run open ipython and test importing from pydrive:

``` python ```
``` from pydrive.auth import GoogleAuth ```

If no error message is thrown. Then you should be good to go.

On your local machine run the python script "localAuthentication.py", to create an authentication key to be used on the server.

``` python localAuthentication.py ```

Type your name when prompted. After Google Authentication, you should see your credentials file in your working directory. Now push these changes to the git.

Now in a separate Terminal window, ssh tunnel into the cans server (try node10), pull the KosikLab repo. Change directory to this folder "pyDrive". Run a test script to be sure that you can access your Google Drive.

Run the testAuthentication.py script:

``` python testAuthentication.py ```

Wait a minute. Should see a "Hello World" txt document appear in your Google Drive. You are all set. Check documentation for syntax when working with Google Drive API [here](https://pythonhosted.org/PyDrive/filemanagement.html)

Be sure to delete your credentials after, if you don't want people to access your Google Drive.