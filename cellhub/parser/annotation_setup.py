import os


class AnnotationSetup:
    def __init__(self, config: dict):
        workdir = config.get("workdir", None)

        # Main annotation directory
        if workdir is not None:
            self.out_dir = os.path.join(workdir, "annotation.dir")
        else:
            self.out_dir = "annotation.dir"

        # API directory for soft links
        if workdir is not None:
            self.api_dir = os.path.join(workdir, "api", "annotation")
        else:
            self.api_dir = os.path.join("api", "annotation")

        # Resource allocations
        default_resources = {
            "threads": 1,
            "mem_mb": 8000,
            "time": "00:10:00",
            "partition": "short",
        }
        self.resources = default_resources.copy()
        self.resources.update(config.get("resources", {}))

        # Annotation parameters
        self.species = config["annotation"]["species"]
        self.ensembl_release = config["annotation"]["ensembl_release"]
        ensembl_host = config["annotation"].get("ensembl_host", "default")
        if ensembl_host == "default":
            self.ensembl_host = ""
        else:
            self.ensembl_host = f"--ensemblhost={ensembl_host}"
