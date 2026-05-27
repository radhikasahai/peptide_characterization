(function () {
  "use strict";

  var API = window.PEPTIDE_API_URL;
  var debounceTimer = null;

  var GEN_MODE_HINTS = {
    random:
      "Builds a sequence of the given length using the standard 20 amino acids (uniform random).",
    motif:
      'Repeats a short motif with an optional N-terminal flank prepended once. Example: motif "RGD", repeats 2, flank "A" → ARGDRGD.',
    combinatorial:
      "Starts from a base sequence and swaps one position (1-based) through every amino acid in Allowed AAs at that position (up to 50 variants). The first variant is loaded into the sequence box.",
  };

  var els = {
    introBenchmark: document.getElementById("intro-benchmark"),
    status: document.getElementById("status"),
    apiConfig: document.getElementById("api-config"),
    apiChangeLink: document.getElementById("api-change-link"),
    apiChangeBtn: document.getElementById("api-change-btn"),
    apiUrlInput: document.getElementById("api-url-input"),
    apiSaveBtn: document.getElementById("api-save-btn"),
    benchmarkSelect: document.getElementById("benchmark-select"),
    sequenceInput: document.getElementById("sequence-input"),
    hydrolyze: document.getElementById("hydrolyze"),
    generateBtn: document.getElementById("generate-btn"),
    genRandom: document.getElementById("gen-random"),
    genMotif: document.getElementById("gen-motif"),
    genCombinatorial: document.getElementById("gen-combinatorial"),
    genModeHint: document.getElementById("gen-mode-hint"),
    libraryCaption: document.getElementById("library-caption"),
    libraryCaptionHint: document.getElementById("library-caption-hint"),
    librarySummary: document.getElementById("library-summary"),
    composition: document.getElementById("composition"),
    compositionSummary: document.getElementById("composition-summary"),
    smiles: document.getElementById("smiles"),
    structureImg: document.getElementById("structure-img"),
    descriptorsBody: document.getElementById("descriptors-body"),
    apiUrl: document.getElementById("api-url"),
  };

  function apiIsConfigured() {
    return API && API.indexOf("__PEPTIDE_API_URL__") === -1;
  }

  function refreshApiDisplay() {
    els.apiUrl.textContent = apiIsConfigured() ? API : "(not set)";
    if (apiIsConfigured()) {
      els.apiChangeLink.classList.remove("hidden");
    } else {
      els.apiChangeLink.classList.add("hidden");
    }
  }

  refreshApiDisplay();

  function formatApiError(data) {
    var detail = data && data.detail;
    if (!detail) {
      return data && data.message ? data.message : "Request failed";
    }
    if (typeof detail === "string") {
      return detail;
    }
    var msg = detail.message || "Request failed";
    if (detail.invalid_residues && detail.invalid_residues.length) {
      msg += " Unsupported letters: " + detail.invalid_residues.join(", ") + ".";
    }
    if (detail.errors && detail.errors.length) {
      msg += " " + detail.errors.join("; ");
    }
    return msg;
  }

  function showStatus(message, type) {
    els.status.textContent = message;
    els.status.className = "status " + (type || "");
    els.status.classList.remove("hidden");
  }

  function hideStatus() {
    els.status.classList.add("hidden");
  }

  function showApiConfig() {
    els.apiConfig.classList.remove("hidden");
    if (apiIsConfigured()) {
      els.apiUrlInput.value = API;
    }
  }

  function hideApiConfig() {
    els.apiConfig.classList.add("hidden");
  }

  function setApiUrl(url) {
    API = url.replace(/\/$/, "");
    window.PEPTIDE_API_URL = API;
    try {
      localStorage.setItem("PEPTIDE_API_URL", API);
    } catch (e) {}
    refreshApiDisplay();
    hideApiConfig();
  }

  function testApiConnection() {
    return apiFetch("/health").then(function (data) {
      if (data.status !== "ok") {
        throw new Error("API health check failed");
      }
    });
  }

  function apiFetch(path, options) {
    if (!apiIsConfigured()) {
      return Promise.reject(
        new Error("Set your Render API URL in the Connect API panel above.")
      );
    }
    return fetch(API + path, options)
      .then(function (response) {
        return response.json().then(function (data) {
          if (!response.ok) {
            throw new Error(formatApiError(data));
          }
          return data;
        });
      })
      .catch(function (err) {
        if (err instanceof TypeError) {
          throw new Error(
            "Cannot reach API at " +
              API +
              ". The service may be waking up—wait ~30s and retry, or check your Render URL."
          );
        }
        throw err;
      });
  }

  function selectedGenMode() {
    var checked = document.querySelector('input[name="gen-mode"]:checked');
    return checked ? checked.value : "random";
  }

  function updateGenPanels() {
    var mode = selectedGenMode();
    els.genRandom.classList.toggle("hidden", mode !== "random");
    els.genMotif.classList.toggle("hidden", mode !== "motif");
    els.genCombinatorial.classList.toggle("hidden", mode !== "combinatorial");
    if (els.genModeHint) {
      els.genModeHint.textContent = GEN_MODE_HINTS[mode] || "";
    }
  }

  function renderCompositionSummary(composition) {
    if (!composition || !els.compositionSummary) {
      return;
    }
    var parts = [
      "Length " + composition.length,
      composition.hydrophobic + " hydrophobic",
      composition.charged + " charged (" +
        composition.positive +
        "+, " +
        composition.negative +
        "−)",
      composition.polar + " polar",
      composition.aromatic + " aromatic",
    ];
    els.compositionSummary.textContent = parts.join(" · ");
    els.compositionSummary.classList.remove("hidden");
  }

  function renderComposition(composition) {
    renderCompositionSummary(composition);
    els.composition.textContent = JSON.stringify(composition, null, 2);
  }

  function clearComposition() {
    if (els.compositionSummary) {
      els.compositionSummary.classList.add("hidden");
      els.compositionSummary.textContent = "";
    }
  }

  function renderDescriptors(descriptors) {
    els.descriptorsBody.innerHTML = "";
    if (!descriptors) {
      var failRow = document.createElement("tr");
      var failCell = document.createElement("td");
      failCell.colSpan = 2;
      failCell.className = "empty-cell";
      failCell.textContent = "Descriptor calculation failed";
      failRow.appendChild(failCell);
      els.descriptorsBody.appendChild(failRow);
      return;
    }
    Object.keys(descriptors)
      .sort()
      .forEach(function (key) {
        var row = document.createElement("tr");
        var keyCell = document.createElement("td");
        var valCell = document.createElement("td");
        keyCell.textContent = key;
        var val = descriptors[key];
        if (typeof val === "number") {
          val = Number.isInteger(val) ? String(val) : val.toFixed(4);
        } else {
          val = val == null ? "" : String(val);
        }
        valCell.textContent = val;
        row.appendChild(keyCell);
        row.appendChild(valCell);
        els.descriptorsBody.appendChild(row);
      });
  }

  function renderResult(data) {
    renderComposition(data.composition);
    els.smiles.textContent = data.smiles || "—";
    if (data.structure_png_base64) {
      els.structureImg.src = "data:image/png;base64," + data.structure_png_base64;
      els.structureImg.classList.remove("hidden");
    } else {
      els.structureImg.classList.add("hidden");
    }
    renderDescriptors(data.descriptors);
    hideStatus();
  }

  function characterize() {
    var sequence = els.sequenceInput.value.trim();
    if (!sequence) {
      clearComposition();
      showStatus("Enter a sequence.", "error");
      return;
    }

    showStatus("Building peptide…", "loading");

    apiFetch("/api/characterize", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        sequence: sequence,
        hydrolyze_sidechains: els.hydrolyze.checked,
      }),
    })
      .then(renderResult)
      .catch(function (err) {
        clearComposition();
        showStatus(err.message, "error");
      });
  }

  function scheduleCharacterize() {
    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(characterize, 350);
  }

  function loadBenchmark() {
    apiFetch("/api/benchmark")
      .then(function (data) {
        var summary = data.summary;
        if (els.introBenchmark) {
          els.introBenchmark.textContent =
            "Regression benchmark loaded: " +
            summary.valid +
            " valid and " +
            summary.invalid +
            " invalid sequences in data/benchmark_sequences.csv.";
          els.introBenchmark.classList.remove("hidden");
        }

        data.examples.forEach(function (ex) {
          var opt = document.createElement("option");
          opt.value = ex.sequence;
          opt.textContent = ex.id + " — " + (ex.notes || ex.sequence);
          opt.dataset.id = ex.id;
          els.benchmarkSelect.appendChild(opt);
        });
      })
      .catch(function (err) {
        showStatus("Could not load benchmark examples: " + err.message, "error");
      });
  }

  function generateSequence() {
    var mode = selectedGenMode();
    var payload = {
      mode: mode,
      hydrolyze_sidechains: els.hydrolyze.checked,
      include_library_summary: mode === "combinatorial",
    };

    if (mode === "random") {
      payload.length = parseInt(document.getElementById("gen-length").value, 10);
    } else if (mode === "motif") {
      payload.motif = document.getElementById("gen-motif-text").value;
      payload.repeats = parseInt(document.getElementById("gen-repeats").value, 10);
      payload.n_term_flank = document.getElementById("gen-flank").value;
    } else {
      payload.base = document.getElementById("gen-base").value;
      payload.position = parseInt(document.getElementById("gen-position").value, 10);
      payload.allowed = document.getElementById("gen-allowed").value;
    }

    showStatus("Generating sequence…", "loading");
    els.generateBtn.disabled = true;

    apiFetch("/api/generate", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    })
      .then(function (data) {
        els.sequenceInput.value = data.sequence;
        els.benchmarkSelect.value = "";

        if (data.variants && data.variants.length) {
          els.libraryCaption.textContent =
            "Library: " + data.variants.length + " sequences (first loaded in input).";
          els.libraryCaption.classList.remove("hidden");
          els.libraryCaptionHint.classList.remove("hidden");
        } else {
          els.libraryCaption.classList.add("hidden");
          els.libraryCaptionHint.classList.add("hidden");
        }

        if (data.library_summary) {
          els.librarySummary.textContent = JSON.stringify(data.library_summary, null, 2);
          els.librarySummary.classList.remove("hidden");
        } else {
          els.librarySummary.classList.add("hidden");
        }

        characterize();
      })
      .catch(function (err) {
        showStatus(err.message, "error");
      })
      .finally(function () {
        els.generateBtn.disabled = false;
      });
  }

  els.apiSaveBtn.addEventListener("click", function () {
    var url = els.apiUrlInput.value.trim();
    if (!url) {
      showStatus("Enter your Render API URL.", "error");
      return;
    }
    if (!/^https?:\/\//.test(url)) {
      showStatus("URL must start with http:// or https://", "error");
      return;
    }
    showStatus("Connecting to API…", "loading");
    setApiUrl(url);
    testApiConnection()
      .then(function () {
        hideStatus();
        loadBenchmark();
        characterize();
      })
      .catch(function (err) {
        showApiConfig();
        showStatus("Could not reach API: " + err.message, "error");
      });
  });

  if (els.apiChangeBtn) {
    els.apiChangeBtn.addEventListener("click", function () {
      showApiConfig();
      els.apiUrlInput.focus();
    });
  }

  els.benchmarkSelect.addEventListener("change", function () {
    if (els.benchmarkSelect.value) {
      els.sequenceInput.value = els.benchmarkSelect.value;
      scheduleCharacterize();
    }
  });

  els.sequenceInput.addEventListener("input", scheduleCharacterize);
  els.hydrolyze.addEventListener("change", scheduleCharacterize);
  els.generateBtn.addEventListener("click", generateSequence);

  document.querySelectorAll('input[name="gen-mode"]').forEach(function (radio) {
    radio.addEventListener("change", updateGenPanels);
  });

  updateGenPanels();
  if (apiIsConfigured()) {
    hideApiConfig();
    loadBenchmark();
    characterize();
  } else {
    showApiConfig();
    showStatus(
      "Deploy the API on Render (see README), then paste your API URL above.",
      "error"
    );
  }
})();
