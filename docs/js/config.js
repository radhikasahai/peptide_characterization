// Set PEPTIDE_API_URL to your deployed API (e.g. Render) before using GitHub Pages.
// Local dev defaults to http://localhost:8000 when opened via file:// or localhost.
(function () {
  if (window.PEPTIDE_API_URL) return;

  var host = window.location.hostname;
  if (host === "localhost" || host === "127.0.0.1") {
    window.PEPTIDE_API_URL = "http://localhost:8000";
    return;
  }

  // Replace with your Render/Railway/Fly API URL after deploying api/
  window.PEPTIDE_API_URL = "https://peptide-characterization-api.onrender.com";
})();
