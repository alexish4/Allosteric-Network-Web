// Import necessary components
import React from 'react';
import { BrowserRouter as Router, Route, Routes, Link, useLocation } from 'react-router-dom';
import Home from './Home';
import Allosteric from './Allosteric'; 

function App() {
  return (
    <Router>
      <Routes>
        <Route path="/" element={<Home />} />
        
        <Route path="/allosteric" element={<Allosteric />} />
      </Routes>
      <NavBar />
    </Router>
  );
}

function NavBar() {
  const location = useLocation(); // Get the current URL

  return (
    <nav>
      {location.pathname !== "/" && <Link to="/">Go Back to Home Page</Link>}
      {location.pathname !== "/allosteric" && <Link to="/allosteric">Allosteric Page</Link>}
    </nav>
  );
}

export default App;