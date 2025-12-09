document.addEventListener('DOMContentLoaded', () => {
    const downloadLink = document.getElementById('downloadLink');
    const messageElement = document.getElementById('message');

    // Add an event listener to the link
    downloadLink.addEventListener('click', () => {
        // Display a confirmation message after the click
        messageElement.textContent = 'Download started! Thank you for downloading.';
        
        // Optionally, clear the message after a few seconds
        setTimeout(() => {
            messageElement.textContent = '';
        }, 5000); 
    });
});