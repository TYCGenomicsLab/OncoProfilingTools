$(document).on("click", ".upload-dropzone, .file-upload-zone, .dataset-upload-zone, .upload-area", function(event) {
  if ($(event.target).is("input[type='file']")) return;

  const fileInput = document.querySelector("#dataset");
  if (fileInput) {
    fileInput.click();
  }
});

$(document).on("click", "label[for='dataset']", function(event) {
  event.preventDefault();

  const fileInput = document.querySelector("#dataset");
  if (fileInput) {
    fileInput.click();
  }
});
