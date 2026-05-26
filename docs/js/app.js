(function () {
  "use strict";

  var API = window.PEPTIDE_API_URL;
  var debounceTimer = null;

  if (!API || API.indexOf("__PEPTIDE_API_URL__") !== -1) {
    document.addEventListener("DOMContentLoaded", function () {
      showStatus(
        "API not configured. Deploy the backend on Render, set the PEPTIDE_API_URL repository variable, and re-run the Pages deploy workflow.",
        "error"
      );
    });
  }

  var els = {
    intro: document.getElementById("intro"),
    status: document.getElementById("status"),
    benchmarkSelect: document.getElementById("benchmark-select"),
    sequenceInput: document.getElementById("sequence-input"),
    hydrolyze: document.getElementById("hydrolyze"),
    generateBtn: document.getElementById("generate-btn"),
    genRandom: document.getElementById("gen-random"),
    genMotif: document.getElementById("gen-motif"),
    genCombinatorial: document.getElementById("gen-combinatorial"),
    libraryCaption: document.getElementById("library-caption"),
    librarySummary: document.getElementById("library-summary"),
    composition: document.getElementById("composition"),
    smiles: document.getElementById("smiles"),
    structureImg: document.getElementById("structure-img"),
    descriptorsBody: document.getElementById("descriptors-body"),
    apiUrl: document.getElementById("api-url"),
  };

  els.apiUrl.textContent = API;

  function showStatus(message, type) {
    els.status.textContent = message;
    els.status.className = "status " + (type || "");
    els.status.classList.remove("hidden");
  }

  function hideStatus() {
    els.status.classList.add("hidden");
  }

  function apiFetch(path, options) {
    if (!API || API.indexOf("__PEPTIDE_API_URL__") !== -1) {
      return Promise.reject(
        new Error("API URL not configured. See README → Web UI → Deploy.")
      );
    }
    return fetch(API + path, options)
      .then(function (response) {
        return response.json().then(function (data) {
          if (!response.ok) {
            var msg =
              (data.detail &&
                (data.detail.message || JSON.stringify(data.detail))) ||
              data.message ||
              "Request failed";
            throw new Error(msg);
          }
          return data;
        });
      })
      .catch(function (err) {
        if (err instanceof TypeError) {
          throw new Error(
            "Cannot reach API at " +
              API +
              ". Deploy the backend on Render and set PEPTIDE_API_URL."
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
  }

  function renderComposition(composition) {
    els.composition.textContent = JSON.stringify(composition, null, 2);
  }

  function renderDescriptors(descriptors) {
    els.descriptorsBody.innerHTML = "";
    if (!descriptors) {
      els.descriptorsBody.innerHTML =
        '<tr><td colspan="2" style="color: var(--muted)">Descriptor calculation failed</td></tr>';
      return;
    }
    Object.keys(descriptors)
      .sort()
      .forEach(function (key) {
        var row = document.createElement("tr");
        var val = descriptors[key];
        if (typeof val === "number") {
          val = Number.isInteger(val) ? val : val.toFixed(4);
        }
        row.innerHTML =
          "<td>" + key + "</td><td>" + val + "</td>";
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
        els.intro.textContent =
          "Built with RDKit: validation and composition from sequence, linear peptide " +
          "assembly via amide coupling, optional Asp/Glu side-chain methyl protection " +
          "during assembly. Regression benchmark: " +
          summary.valid +
          " valid and " +
          summary.invalid +
          " invalid sequences.";

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
        } else {
          els.libraryCaption.classList.add("hidden");
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
  loadBenchmark();
  characterize();
})();
