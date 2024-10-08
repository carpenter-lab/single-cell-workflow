#!/usr/bin/env python3

import hashlib
import os
import pathlib
import pickle
import subprocess
import urllib
from typing import Optional, Self

import click
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
                driveId=self._drive_id,
                includeItemsFromAllDrives=True,
                supportsAllDrives=True,
                fields="files(id, name, mimeType, size, parents, modifiedTime)",
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
            file_id = self.get_file_id(filename)
            self.validate_service()
            self.service.files().update(
                supportsAllDrives=True,
                fileId=file_id,
                media_body=media,
                body=metadata,
            ).execute()
        except HttpError as e:
            raise RuntimeWarning(
                f"Failed to update file: {filename} because {e.resp['message']}"
            )


@click.group()
def cli():
    pass


@cli.command()
def dry_run():
    main(["-c", "1", "--dry-run", "--quiet", "rules"])


@cli.command()
def setup_google_drive():
    if not os.path.exists("client_secrets.json"):
        raise RuntimeError("A `client_secrets.json` file is required.")
    UploadToGoogleDrive()


@cli.command()
@click.option("--cores", default=1, help="Number of cores", type=click.INT)
@click.option("--report", default="report.html", help="Path to use to generate report", type=click.Path(dir_okay=False, path_type=pathlib.Path))
@click.option("--gdrive-name", default="report.html", help="Name to use for report on Google Drive", type=click.STRING)
@click.option("--gdrive", default=None, help="Drive to use in Google Drive", type=click.STRING)
def run_workflow(cores, report, gdrive_name, gdrive):
    try:
        import panoptes
        urllib.request.urlopen("http://127.0.0.1:5000", timeout=1)
        opts = ["--wms-monitor", "http://127.0.0.1:5000"]
    except ModuleNotFoundError:
        opts = []
    except urllib.error.URLError:
        opts = []

    try:
        sbatch_ready = subprocess.run(["sbatch", "--version"], stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        sbatch_ready = False

    if pathlib.Path("profiles/slurm").exists() and sbatch_ready:
        profile_args = ["--profile", "profiles/slurm"]
    elif pathlib.Path("profiles/local").exists():
        profile_args = ["--profile", "profiles/local"]
    else:
        profile_args = ["-c", str(cores)]

    opts.extend(profile_args)

    try:
        main(opts)
    except SystemExit as e:
        if e.code != 0:
            raise e

    try:
        main(["--report", str(report), "--quiet", "all"])
    except SystemExit:
        pass

    if os.path.exists("client_secrets.json") and gdrive:
        print("Uploading to Google Drive...")
        UploadToGoogleDrive().add_drive_label(gdrive).update_with_new_file(
            gdrive_name, report
        )


if __name__ == "__main__":
    cli(auto_envvar_prefix="WORKFLOW")
