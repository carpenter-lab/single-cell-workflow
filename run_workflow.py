import hashlib
import os
import pathlib
import pickle
from typing import Optional, Self

from google.auth.transport.requests import Request
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import Resource, build
from googleapiclient.errors import HttpError
from googleapiclient.http import MediaFileUpload
from oauthlib.oauth2 import AccessDeniedError
from snakemake.cli import main


class ResourceError(Exception):
    pass


class UploadToGoogleDrive:
    _scopes = ["https://www.googleapis.com/auth/drive"]

    def __init__(self) -> None:
        self._drive_id = None
        self._drive_label = None
        self._service = self.get_gdrive_service()

    @property
    def service(self) -> Resource:
        return self._service

    @property
    def drive_label(self) -> str:
        return self._drive_label

    @drive_label.setter
    def drive_label(self, drive_label) -> None:
        self._drive_label = drive_label
        self.get_drive_id(drive_label)

    def add_drive_label(self, drive_label) -> Self:
        self.drive_label = drive_label
        return self

    def get_gdrive_service(self) -> Optional[Resource]:
        creds = None
        # The file token.pickle stores the user's access and refresh tokens and is
        # created automatically when the authorization flow completes for the first
        # time.
        if os.path.exists("token.pickle"):
            with open("token.pickle", "rb") as token:
                creds = pickle.load(token)
        # If there are no (valid) credentials available, let the user log in.
        if not creds or not creds.valid:
            if creds and creds.expired and creds.refresh_token:
                creds.refresh(Request())
            else:
                flow = InstalledAppFlow.from_client_secrets_file(
                    "client_secrets.json", self._scopes
                )
                creds = flow.run_local_server(port=0)
            # Save the credentials for the next run
            with open("token.pickle", "wb") as token:
                pickle.dump(creds, token)
        # return Google Drive API service

        try:
            return build("drive", "v3", credentials=creds)
        except AccessDeniedError:
            raise ResourceWarning("Failed to Authenticate to Google Drive")

    @classmethod
    def validate_service(cls) -> None:
        if not cls.service:
            raise ResourceError("Service is missing")

    def get_drive_id(self, drive_label) -> None:
        self.validate_service()
        drives = self.service.drives().list().execute()

        for drive in drives["drives"]:
            if drive["name"] == drive_label:  # Robert Drive ID: 0APjpJCIfwQE2Uk9PVA
                self._drive_id = drive["id"]

    def get_file_id(self, filename) -> Optional[str]:
        """Function to get the file id for a given filename.

        :param filename: the name of the file to search for.
        :return: the file id if found, None otherwise.

        """
        self.validate_service()
        files = (
            self.service.files()
            .list(
                corpora="drive",
                pageSize=5,
                driveId=self._drive_id,
                includeItemsFromAllDrives=True,
                supportsAllDrives=True,
                fields="nextPageToken, files(id, name, mimeType, size, parents, modifiedTime)",
            )
            .execute()
        )

        for file in files["files"]:
            if file["name"] == filename:
                return file["id"]

    def update_with_new_file(self, filename: str, local_file: str) -> None:
        """
        Updates an existing file with the contents of a local file.

        :param filename: The name of the file to be updated.
        :type filename: str
        :param local_file: The path of the local file to upload.
        :type local_file: str
        :return: None
        :rtype: None
        """

        media = MediaFileUpload(local_file)
        metadata = {"name": filename, "mimetype": "text/html"}

        try:
            self.validate_service()
            self.service.files().update(
                supportsAllDrives=True,
                fileId=self.get_file_id(filename),
                media_body=media,
                body=metadata,
            ).execute()
        except HttpError as e:
            raise RuntimeWarning(
                f"Failed to update file: {filename} because {e.resp['message']}"
            )


if __name__ == "__main__":
    if os.getenv("DRYRUN") == "1":
        main(["-c", "1", "--dry-run", "--quiet", "rules"])
    if os.getenv("SETUP_GDRIVE") == 1:
        if not os.path.exists("client_secrets.json"):
            raise RuntimeError("A `client_secrets.json` file is required.")
        UploadToGoogleDrive()
    else:
        try:
            main(["-c", os.getenv("WORKFLOW_CORES")])
        except SystemExit as e:
            if e.code != 0:
                raise e

        hash_old = None
        if pathlib.Path("report.txt").exists():
            with open("report.html", "rb") as f:
                hash_old = hashlib.md5(f.read()).hexdigest()

        try:
            main(["--report", "report.html", "--quiet", "all"])
        except SystemExit as e:
            if e.code != 0:
                raise e

        with open("report.html", "rb") as f:
            hash_new = hashlib.md5(f.read()).hexdigest()

        if (hash_old != hash_new) & os.path.exists("client_secrets.json"):
            UploadToGoogleDrive().add_drive_label("Robert").update_with_new_file(
                "BAL scRNA seq.html", "report.html"
            )
